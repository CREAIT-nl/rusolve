//! The error type of the `fn ... -> Result<_, CoreError>` cores that the CLI
//! subcommands (`explain`, `lint`, `ranging`, `feasrelax`) and the MCP server
//! (`mcp.rs`) share. The CLI turns it into an exit code (exit-code table in
//! `mip_highs.rs`) and a stderr line; the MCP server turns it into a tool result
//! with `isError: true`.

/// A failure of a core: the exit code the CLI exits with, and the message
/// the CLI prints after its `rusolve <command>:` prefix.
#[derive(Debug)]
pub struct CoreError {
    /// 2 usage, 3 input, 4 internal.
    pub code: u8,
    pub msg: String,
}

impl CoreError {
    pub fn usage(msg: impl Into<String>) -> Self {
        Self {
            code: 2,
            msg: msg.into(),
        }
    }
    pub fn internal(msg: impl Into<String>) -> Self {
        Self {
            code: 4,
            msg: msg.into(),
        }
    }
}
