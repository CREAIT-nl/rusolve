//! The skills as MCP prompts, the schemas and the
//! agent reference as MCP resources.
//!
//! Prompts: <https://modelcontextprotocol.io/specification/2026-07-28/server/prompts>
//! (`prompts/list`, `prompts/get`: a prompt is a name, a description and the messages it expands
//! to). Resources: <https://modelcontextprotocol.io/specification/2026-07-28/server/resources>
//! (`resources/list`, `resources/read`; an unknown URI is error -32002).
//!
//! Everything is embedded at build time from `docs/skills/*/SKILL.md`, `docs/schema/*.json` and
//! `docs/agent/reference.md` (`build.rs` writes `$OUT_DIR/mcp_docs.rs`), so a skill or schema added
//! to the tree appears in the next build with no Rust change. A prompt's name is its directory
//! name (`solve`, `debug-infeasible`, ...); its `title` is the skill's own `name`, its
//! `description` the skill's `description`, and the message is the skill body without the YAML
//! front matter. A skill takes no arguments.

use serde_json::{Value, json};

include!(concat!(env!("OUT_DIR"), "/mcp_docs.rs"));

/// JSON-RPC code for an unknown resource URI (resources spec, "Error Handling").
pub const RESOURCE_NOT_FOUND: i64 = -32002;

const URI_SCHEMA: &str = "rusolve://schema/";
const URI_REFERENCE: &str = "rusolve://docs/agent/reference";

/// `(front matter value of `key`, body)` of a SKILL.md: the `---` block at the top holds
/// `key: value` lines; the body is what follows it.
fn split_front_matter(text: &str) -> (Vec<(String, String)>, &str) {
    let Some(rest) = text.strip_prefix("---\n") else {
        return (Vec::new(), text);
    };
    let Some(end) = rest.find("\n---") else {
        return (Vec::new(), text);
    };
    let fields = rest[..end]
        .lines()
        .filter_map(|l| {
            let (k, v) = l.split_once(':')?;
            Some((k.trim().to_string(), v.trim().to_string()))
        })
        .collect();
    let body = rest[end + 4..].trim_start_matches('\n');
    (fields, body)
}

fn field<'a>(fields: &'a [(String, String)], key: &str) -> Option<&'a str> {
    fields
        .iter()
        .find(|(k, _)| k == key)
        .map(|(_, v)| v.as_str())
}

pub fn prompts_list() -> Value {
    let prompts: Vec<Value> = SKILLS
        .iter()
        .map(|(name, text)| {
            let (fm, _) = split_front_matter(text);
            json!({
                "name": name,
                "title": field(&fm, "name").unwrap_or(name),
                "description": field(&fm, "description").unwrap_or(""),
                "arguments": [],
            })
        })
        .collect();
    json!({"prompts": prompts})
}

/// `prompts/get`: `Err` is the message of a -32602 (unknown prompt name).
pub fn prompts_get(name: &str) -> Result<Value, String> {
    let Some((_, text)) = SKILLS.iter().find(|(n, _)| *n == name) else {
        let known: Vec<&str> = SKILLS.iter().map(|(n, _)| *n).collect();
        return Err(format!(
            "Unknown prompt: {name} (prompts/list names: {})",
            known.join(", ")
        ));
    };
    let (fm, body) = split_front_matter(text);
    Ok(json!({
        "description": field(&fm, "description").unwrap_or(""),
        "messages": [{"role": "user", "content": {"type": "text", "text": body}}],
    }))
}

fn schema_stem(file: &str) -> &str {
    file.strip_suffix(".schema.json")
        .or_else(|| file.strip_suffix(".json"))
        .unwrap_or(file)
}

pub fn resources_list() -> Value {
    let mut resources: Vec<Value> = SCHEMAS
        .iter()
        .map(|(file, text)| {
            let stem = schema_stem(file);
            let desc = serde_json::from_str::<Value>(text)
                .ok()
                .and_then(|v| v.get("title").and_then(Value::as_str).map(str::to_string))
                .unwrap_or_else(|| format!("JSON Schema of the {stem} object"));
            json!({
                "uri": format!("{URI_SCHEMA}{stem}"),
                "name": stem,
                "title": format!("Schema {stem}"),
                "description": desc,
                "mimeType": "application/schema+json",
                "size": text.len(),
            })
        })
        .collect();
    resources.push(json!({
        "uri": URI_REFERENCE,
        "name": "agent-reference",
        "title": "rusolve agent reference",
        "description": "Every command, flag and option of rusolve, generated from the registries (docs/agent/reference.md).",
        "mimeType": "text/markdown",
        "size": AGENT_REFERENCE.len(),
    }));
    json!({"resources": resources})
}

/// `resources/read`: `Err` is the message of a -32002 (unknown URI).
pub fn resources_read(uri: &str) -> Result<Value, String> {
    let (mime, text) = if uri == URI_REFERENCE {
        ("text/markdown", AGENT_REFERENCE)
    } else if let Some(stem) = uri.strip_prefix(URI_SCHEMA) {
        match SCHEMAS.iter().find(|(f, _)| schema_stem(f) == stem) {
            Some((_, t)) => ("application/schema+json", *t),
            None => return Err(format!("Resource not found: {uri}")),
        }
    } else {
        return Err(format!("Resource not found: {uri}"));
    };
    Ok(json!({"contents": [{"uri": uri, "mimeType": mime, "text": text}]}))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn front_matter_splits() {
        let (fm, body) = split_front_matter("---\nname: a\ndescription: b: c\n---\n\n# H\n");
        assert_eq!(field(&fm, "name"), Some("a"));
        assert_eq!(field(&fm, "description"), Some("b: c"));
        assert_eq!(body, "# H\n");
    }

    #[test]
    fn every_embedded_skill_is_a_prompt() {
        let list = prompts_list();
        assert_eq!(list["prompts"].as_array().unwrap().len(), SKILLS.len());
        assert!(!SKILLS.is_empty());
        for (n, _) in SKILLS {
            assert!(prompts_get(n).is_ok());
        }
    }
}
