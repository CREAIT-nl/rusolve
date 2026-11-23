//! HiGHS source: util/HighsRbTree.h (lines 1-452)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.
//!
//! C++ uses CRTP: `class X : public RbTree<X>` / `CacheMinRbTree<X>`, where
//! `X` supplies `getRbTreeLinks(node)` and `getKey(node)` and the traits
//! `RbTreeTraits<X>::{KeyType, LinkType}`; the tree holds `LinkType&
//! rootNode` (and `LinkType& first_`) into the owner. In Rust the CRTP base
//! is the trait [`RbTree`] (provided methods = the C++ `RbTree<Impl>`
//! members) and [`CacheMinRbTree`] (the cached-minimum overrides of `link`/
//! `unlink`/`first`). The implementing types are *view structs* holding
//! `&'a mut` borrows of the owner's node storage, root and first fields
//! (see `HighsNodeQueue_NodeLowerRbTree`), constructed per use exactly where
//! C++ constructs its temporary tree object (`NodeLowerRbTree rbTree(this);`).
//!
//! `KeyType` comparisons use only `<` (`PartialOrd::lt`), as C++ does;
//! tuple `<` is lexicographic like `std::tuple::operator<` for non-NaN keys.

use crate::prelude::*;

/// Link types used by HiGHS red-black trees (`HighsInt` and `int64_t`).
/// `NO_LINK` = `RbTreeLinks<T>::noLink()` (= -1); `STORAGE_BITS` =
/// `sizeof(ParentStorageType) * 8` (`make_unsigned<T>`).
pub trait RbLink: Copy + PartialEq + std::fmt::Debug + Default {
    const NO_LINK: Self;
    const STORAGE_BITS: u32;
    fn to_storage(self) -> u64;
    fn from_storage(v: u64) -> Self;
}

impl RbLink for i32 {
    const NO_LINK: Self = -1;
    const STORAGE_BITS: u32 = 32;
    fn to_storage(self) -> u64 {
        // (ParentStorageType)(value) for the unsigned 32-bit storage type
        u64::from(self as u32)
    }
    fn from_storage(v: u64) -> Self {
        // (T)(storage) truncates to the 32-bit type
        v as u32 as i32
    }
}

impl RbLink for i64 {
    const NO_LINK: Self = -1;
    const STORAGE_BITS: u32 = 64;
    fn to_storage(self) -> u64 {
        self as u64
    }
    fn from_storage(v: u64) -> Self {
        v as i64
    }
}

/// HiGHS: HighsRbTree.h::highs::RbTreeLinks<T> (lines 19-84)
///
/// C++: `T child[2]; ParentStorageType parentAndColor;` The parent is stored
/// as `p + 1` with the colour in the top bit of the unsigned storage type
/// (`colorBitPos() = sizeof(ParentStorageType) * 8 - 1`); kept in a `u64`
/// here, using only the low `T::STORAGE_BITS` bits, so the encoding is the same.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RbTreeLinks<T: RbLink> {
    /// C++: `T child[2]` (`kLeft = 0`, `kRight = 1`)
    pub child: [T; 2],
    /// C++: `ParentStorageType parentAndColor`
    pub parent_and_color: u64,
}

impl<T: RbLink> RbTreeLinks<T> {
    pub const kLeft: usize = 0;
    pub const kRight: usize = 1;

    /// HiGHS: HighsRbTree.h::RbTreeLinks::colorBitPos (lines 32-34)
    pub fn color_bit_pos() -> u32 {
        T::STORAGE_BITS - 1
    }
    /// HiGHS: HighsRbTree.h::RbTreeLinks::colorBitMask (lines 36-38)
    pub fn color_bit_mask() -> u64 {
        1u64 << Self::color_bit_pos()
    }
    /// HiGHS: HighsRbTree.h::RbTreeLinks::noLink (lines 42-52)
    pub fn no_link() -> T {
        T::NO_LINK
    }
    /// HiGHS: HighsRbTree.h::RbTreeLinks::getColor (line 56)
    pub fn get_color(&self) -> bool {
        (self.parent_and_color >> Self::color_bit_pos()) != 0
    }
    /// HiGHS: HighsRbTree.h::RbTreeLinks::isBlack (line 58)
    pub fn is_black(&self) -> bool {
        !self.get_color()
    }
    /// HiGHS: HighsRbTree.h::RbTreeLinks::isRed (line 60)
    pub fn is_red(&self) -> bool {
        self.get_color()
    }
    /// HiGHS: HighsRbTree.h::RbTreeLinks::makeRed (line 62)
    pub fn make_red(&mut self) {
        self.parent_and_color |= Self::color_bit_mask();
    }
    /// HiGHS: HighsRbTree.h::RbTreeLinks::makeBlack (line 64)
    pub fn make_black(&mut self) {
        self.parent_and_color &= !Self::color_bit_mask();
    }
    /// HiGHS: HighsRbTree.h::RbTreeLinks::setColor (lines 64-67)
    pub fn set_color(&mut self, color: bool) {
        self.make_black();
        self.parent_and_color |= u64::from(color) << Self::color_bit_pos();
    }
    /// HiGHS: HighsRbTree.h::RbTreeLinks::setParent (lines 69-76)
    pub fn set_parent(&mut self, p: T) {
        // the storage holds only the low STORAGE_BITS bits (ParentStorageType)
        let storage_mask: u64 = if T::STORAGE_BITS >= 64 {
            u64::MAX
        } else {
            (1u64 << T::STORAGE_BITS) - 1
        };
        // (ParentStorageType)(p + 1)
        self.parent_and_color = (Self::color_bit_mask() & self.parent_and_color)
            | (p.to_storage().wrapping_add(1) & storage_mask);
    }
    /// HiGHS: HighsRbTree.h::RbTreeLinks::getParent (lines 78-83)
    pub fn get_parent(&self) -> T {
        // ((T)(parentAndColor & ~colorBitMask())) - 1
        T::from_storage((self.parent_and_color & !Self::color_bit_mask()).wrapping_sub(1))
    }
}

/// HiGHS: HighsRbTree.h::highs::RbTree<Impl> (lines 89-418), CRTP base as a trait.
///
/// Required methods = what `Impl` and `RbTreeTraits<Impl>` supply in C++
/// (`getRbTreeLinks`, `getKey`, and the `LinkType& rootNode` reference).
/// Provided methods = the C++ `RbTree` members.
/// C++ `Dir` (`kLeft = 0`, `kRight = 1`) is `usize` here.
pub trait RbTree {
    /// C++: `RbTreeTraits<Impl>::KeyType`
    type KeyType: PartialOrd;
    /// C++: `RbTreeTraits<Impl>::LinkType`
    type LinkType: RbLink;

    /// C++: reading `LinkType& rootNode`
    fn root_node(&self) -> Self::LinkType;
    /// C++: writing `LinkType& rootNode`
    fn set_root_node(&mut self, root: Self::LinkType);
    /// C++: `Impl::getRbTreeLinks(node) const`
    fn get_rb_tree_links(&self, node: Self::LinkType) -> &RbTreeLinks<Self::LinkType>;
    /// C++: `Impl::getRbTreeLinks(node)` (non-const)
    fn get_rb_tree_links_mut(&mut self, node: Self::LinkType) -> &mut RbTreeLinks<Self::LinkType>;
    /// C++: `Impl::getKey(node) const`
    fn get_key(&self, node: Self::LinkType) -> Self::KeyType;

    /// HiGHS: HighsRbTree.h::RbTree::isRed (lines 122-125)
    fn is_red(&self, node: Self::LinkType) -> bool {
        node != Self::LinkType::NO_LINK && self.get_rb_tree_links(node).is_red()
    }
    /// HiGHS: HighsRbTree.h::RbTree::isBlack (lines 127-130)
    fn is_black(&self, node: Self::LinkType) -> bool {
        node == Self::LinkType::NO_LINK || self.get_rb_tree_links(node).is_black()
    }
    /// HiGHS: HighsRbTree.h::RbTree::makeRed (lines 132-134)
    fn make_red(&mut self, node: Self::LinkType) {
        self.get_rb_tree_links_mut(node).make_red();
    }
    /// HiGHS: HighsRbTree.h::RbTree::makeBlack (lines 136-138)
    fn make_black(&mut self, node: Self::LinkType) {
        self.get_rb_tree_links_mut(node).make_black();
    }
    /// HiGHS: HighsRbTree.h::RbTree::setColor (lines 140-142)
    fn set_color(&mut self, node: Self::LinkType, color: HighsUInt) {
        self.get_rb_tree_links_mut(node).set_color(color != 0);
    }
    /// HiGHS: HighsRbTree.h::RbTree::getColor (lines 144-146)
    fn get_color(&self, node: Self::LinkType) -> HighsUInt {
        HighsUInt::from(self.get_rb_tree_links(node).get_color())
    }
    /// HiGHS: HighsRbTree.h::RbTree::setParent (lines 148-150)
    fn set_parent(&mut self, node: Self::LinkType, parent: Self::LinkType) {
        self.get_rb_tree_links_mut(node).set_parent(parent);
    }
    /// HiGHS: HighsRbTree.h::RbTree::getParent (lines 152-154)
    fn get_parent(&self, node: Self::LinkType) -> Self::LinkType {
        self.get_rb_tree_links(node).get_parent()
    }
    /// HiGHS: HighsRbTree.h::RbTree::opposite (line 160)
    fn opposite(dir: usize) -> usize
    where
        Self: Sized,
    {
        1 - dir
    }
    /// HiGHS: HighsRbTree.h::RbTree::getChild (lines 162-164)
    fn get_child(&self, node: Self::LinkType, dir: usize) -> Self::LinkType {
        self.get_rb_tree_links(node).child[dir]
    }
    /// HiGHS: HighsRbTree.h::RbTree::setChild (lines 166-168)
    fn set_child(&mut self, node: Self::LinkType, dir: usize, child: Self::LinkType) {
        self.get_rb_tree_links_mut(node).child[dir] = child;
    }
    /// HiGHS: HighsRbTree.h::RbTree::rotate (lines 170-186)
    fn rotate(&mut self, x: Self::LinkType, dir: usize) {
        let no_link = Self::LinkType::NO_LINK;
        // opposite(dir) is written as `1 - dir`: `Self::opposite` needs `Self: Sized`,
        // which a provided method cannot assume.
        let y = self.get_child(x, 1 - dir);
        let y_dir = self.get_child(y, dir);
        self.set_child(x, 1 - dir, y_dir);
        if y_dir != no_link {
            self.set_parent(y_dir, x);
        }

        let p_x = self.get_parent(x);
        self.set_parent(y, p_x);

        if p_x == no_link {
            self.set_root_node(y);
        } else {
            let side = usize::from(x != self.get_child(p_x, dir)) ^ dir;
            self.set_child(p_x, side, y);
        }

        self.set_child(y, dir, x);
        self.set_parent(x, y);
    }
    /// HiGHS: HighsRbTree.h::RbTree::insertFixup (lines 188-220)
    fn insert_fixup(&mut self, z: Self::LinkType) {
        let no_link = Self::LinkType::NO_LINK;
        let k_left = RbTreeLinks::<Self::LinkType>::kLeft;
        let mut z = z;

        let mut p_z = self.get_parent(z);
        while self.is_red(p_z) {
            let mut z_grand_parent = self.get_parent(p_z);
            debug_assert!(z_grand_parent != no_link);

            let dir = usize::from(self.get_child(z_grand_parent, k_left) == p_z);

            let y = self.get_child(z_grand_parent, dir);
            if self.is_red(y) {
                self.make_black(p_z);
                self.make_black(y);
                self.make_red(z_grand_parent);
                z = z_grand_parent;
            } else {
                if z == self.get_child(p_z, dir) {
                    z = p_z;
                    self.rotate(z, 1 - dir);
                    p_z = self.get_parent(z);
                    z_grand_parent = self.get_parent(p_z);
                    debug_assert!(z_grand_parent != no_link);
                }

                self.make_black(p_z);
                self.make_red(z_grand_parent);
                self.rotate(z_grand_parent, dir);
            }

            p_z = self.get_parent(z);
        }

        let root = self.root_node();
        self.make_black(root);
    }
    /// HiGHS: HighsRbTree.h::RbTree::transplant (lines 222-234)
    fn transplant(
        &mut self,
        u: Self::LinkType,
        v: Self::LinkType,
        nil_parent: &mut Self::LinkType,
    ) {
        let no_link = Self::LinkType::NO_LINK;
        let k_left = RbTreeLinks::<Self::LinkType>::kLeft;

        let p = self.get_parent(u);

        if p == no_link {
            self.set_root_node(v);
        } else {
            let side = usize::from(u != self.get_child(p, k_left));
            self.set_child(p, side, v);
        }

        if v == no_link {
            *nil_parent = p;
        } else {
            self.set_parent(v, p);
        }
    }
    /// HiGHS: HighsRbTree.h::RbTree::deleteFixup (lines 236-281)
    fn delete_fixup(&mut self, x: Self::LinkType, nil_parent: Self::LinkType) {
        let no_link = Self::LinkType::NO_LINK;
        let k_left = RbTreeLinks::<Self::LinkType>::kLeft;
        let k_right = RbTreeLinks::<Self::LinkType>::kRight;
        let mut x = x;

        while x != self.root_node() && self.is_black(x) {
            let p = if x == no_link {
                nil_parent
            } else {
                self.get_parent(x)
            };
            let dir = usize::from(x == self.get_child(p, k_left));
            let mut w = self.get_child(p, dir);
            debug_assert!(w != no_link);

            if self.is_red(w) {
                self.make_black(w);
                self.make_red(p);
                self.rotate(p, 1 - dir);
                debug_assert!(
                    (x == no_link && p == nil_parent) || (x != no_link && p == self.get_parent(x))
                );
                w = self.get_child(p, dir);
                debug_assert!(w != no_link);
            }

            if self.is_black(self.get_child(w, k_left)) && self.is_black(self.get_child(w, k_right))
            {
                self.make_red(w);
                x = p;
            } else {
                if self.is_black(self.get_child(w, dir)) {
                    let w_opposite = self.get_child(w, 1 - dir);
                    self.make_black(w_opposite);
                    self.make_red(w);
                    self.rotate(w, dir);
                    debug_assert!(
                        (x == no_link && p == nil_parent)
                            || (x != no_link && p == self.get_parent(x))
                    );
                    w = self.get_child(p, dir);
                }
                let p_color = self.get_color(p);
                self.set_color(w, p_color);
                self.make_black(p);
                let w_dir = self.get_child(w, dir);
                self.make_black(w_dir);
                self.rotate(p, 1 - dir);
                x = self.root_node();
            }
        }

        if x != no_link {
            self.make_black(x);
        }
    }
    /// HiGHS: HighsRbTree.h::RbTree::empty (line 281)
    fn empty(&self) -> bool {
        self.root_node() == Self::LinkType::NO_LINK
    }
    /// HiGHS: HighsRbTree.h::RbTree::first (lines 283-291), `first(LinkType x)`
    fn first_from(&self, x: Self::LinkType) -> Self::LinkType {
        let no_link = Self::LinkType::NO_LINK;
        let k_left = RbTreeLinks::<Self::LinkType>::kLeft;
        let mut x = x;

        if x == no_link {
            return no_link;
        }

        loop {
            let l_x = self.get_child(x, k_left);
            if l_x == no_link {
                return x;
            }
            x = l_x;
        }
    }
    /// HiGHS: HighsRbTree.h::RbTree::first() (line 303, = `first(rootNode)`)
    ///
    /// `rb_` prefix: a [`CacheMinRbTree`] hides this with its cached `first()`
    /// (C++ name hiding). The prefix keeps the base reachable only by an explicit
    /// name, so `tree.first()` on a cache-min tree cannot resolve here.
    fn rb_first(&self) -> Self::LinkType {
        self.first_from(self.root_node())
    }
    /// HiGHS: HighsRbTree.h::RbTree::last (lines 293-301), `last(LinkType x)`
    fn last_from(&self, x: Self::LinkType) -> Self::LinkType {
        let no_link = Self::LinkType::NO_LINK;
        let k_right = RbTreeLinks::<Self::LinkType>::kRight;
        let mut x = x;

        if x == no_link {
            return no_link;
        }

        loop {
            let r_x = self.get_child(x, k_right);
            if r_x == no_link {
                return x;
            }
            x = r_x;
        }
    }
    /// HiGHS: HighsRbTree.h::RbTree::last() (line 305, = `last(rootNode)`)
    fn last(&self) -> Self::LinkType {
        self.last_from(self.root_node())
    }
    /// HiGHS: HighsRbTree.h::RbTree::successor (lines 307-318)
    fn successor(&self, x: Self::LinkType) -> Self::LinkType {
        let no_link = Self::LinkType::NO_LINK;
        let k_right = RbTreeLinks::<Self::LinkType>::kRight;
        let mut x = x;

        let mut y = self.get_child(x, k_right);
        if y != no_link {
            return self.first_from(y);
        }

        y = self.get_parent(x);
        while y != no_link && x == self.get_child(y, k_right) {
            x = y;
            y = self.get_parent(x);
        }

        y
    }
    /// HiGHS: HighsRbTree.h::RbTree::predecessor (lines 320-331)
    fn predecessor(&self, x: Self::LinkType) -> Self::LinkType {
        let no_link = Self::LinkType::NO_LINK;
        let k_left = RbTreeLinks::<Self::LinkType>::kLeft;
        let mut x = x;

        let mut y = self.get_child(x, k_left);
        if y != no_link {
            return self.last_from(y);
        }

        y = self.get_parent(x);
        while y != no_link && x == self.get_child(y, k_left) {
            x = y;
            y = self.get_parent(x);
        }

        y
    }
    /// HiGHS: HighsRbTree.h::RbTree::find (lines 333-352), `find(key, treeRoot)`
    fn find_from(&self, key: &Self::KeyType, tree_root: Self::LinkType) -> (Self::LinkType, bool) {
        let no_link = Self::LinkType::NO_LINK;
        let k_left = RbTreeLinks::<Self::LinkType>::kLeft;
        let k_right = RbTreeLinks::<Self::LinkType>::kRight;

        let mut y = no_link;
        let mut x = tree_root;
        while x != no_link {
            // cmp = 1 - (getKey(x) < key) + (key < getKey(x))
            let x_lt_key = PartialOrd::lt(&self.get_key(x), key);
            let key_lt_x = PartialOrd::lt(key, &self.get_key(x));
            let cmp: HighsInt = 1 - HighsInt::from(x_lt_key) + HighsInt::from(key_lt_x);
            match cmp {
                0 => {
                    y = x;
                    x = self.get_child(y, k_right);
                }
                1 => return (x, true),
                2 => {
                    y = x;
                    x = self.get_child(y, k_left);
                }
                // cmp is always 0, 1 or 2 (C++ switch has no default)
                _ => {}
            }
        }

        (y, false)
    }
    /// HiGHS: HighsRbTree.h::RbTree::find (lines 354-356), `find(key)`
    fn find(&self, key: &Self::KeyType) -> (Self::LinkType, bool) {
        self.find_from(key, self.root_node())
    }
    /// HiGHS: HighsRbTree.h::RbTree::link (lines 358-369), `link(z, parent)`
    /// (`rb_` prefix: see [`RbTree::rb_first`]).
    fn rb_link_parent(&mut self, z: Self::LinkType, parent: Self::LinkType) {
        let no_link = Self::LinkType::NO_LINK;
        let k_left = RbTreeLinks::<Self::LinkType>::kLeft;
        let k_right = RbTreeLinks::<Self::LinkType>::kRight;

        self.set_parent(z, parent);
        if parent == no_link {
            self.set_root_node(z);
        } else {
            let side = usize::from(PartialOrd::lt(&self.get_key(parent), &self.get_key(z)));
            self.set_child(parent, side, z);
        }

        self.set_child(z, k_left, no_link);
        self.set_child(z, k_right, no_link);
        self.make_red(z);
        self.insert_fixup(z);
    }
    /// HiGHS: HighsRbTree.h::RbTree::link (lines 371-380), `link(z)`
    /// for a plain `RbTree` (`Impl::link(z, y)` = [`RbTree::rb_link_parent`]).
    /// A [`CacheMinRbTree`] uses its own `link(z)`, which dispatches to the
    /// cached `link_parent` as C++'s `static_cast<Impl*>(this)->link(z, y)` does.
    fn rb_link(&mut self, z: Self::LinkType) {
        let no_link = Self::LinkType::NO_LINK;

        let mut y = no_link;
        let mut x = self.root_node();
        while x != no_link {
            y = x;
            let side = usize::from(PartialOrd::lt(&self.get_key(x), &self.get_key(z)));
            x = self.get_child(y, side);
        }

        self.rb_link_parent(z, y);
    }
    /// HiGHS: HighsRbTree.h::RbTree::unlink (lines 382-416)
    /// (`rb_` prefix: see [`RbTree::rb_first`]).
    fn rb_unlink(&mut self, z: Self::LinkType) {
        let no_link = Self::LinkType::NO_LINK;
        let k_left = RbTreeLinks::<Self::LinkType>::kLeft;
        let k_right = RbTreeLinks::<Self::LinkType>::kRight;

        let mut nil_parent = no_link;
        let mut y = z;
        let mut y_was_black = self.is_black(y);
        let x: Self::LinkType;

        if self.get_child(z, k_left) == no_link {
            x = self.get_child(z, k_right);
            self.transplant(z, x, &mut nil_parent);
        } else if self.get_child(z, k_right) == no_link {
            x = self.get_child(z, k_left);
            self.transplant(z, x, &mut nil_parent);
        } else {
            y = self.first_from(self.get_child(z, k_right));
            y_was_black = self.is_black(y);
            x = self.get_child(y, k_right);
            if self.get_parent(y) == z {
                if x == no_link {
                    nil_parent = y;
                } else {
                    self.set_parent(x, y);
                }
            } else {
                let y_right = self.get_child(y, k_right);
                self.transplant(y, y_right, &mut nil_parent);
                let z_right = self.get_child(z, k_right);
                self.set_child(y, k_right, z_right);
                self.set_parent(z_right, y);
            }
            self.transplant(z, y, &mut nil_parent);
            let z_left = self.get_child(z, k_left);
            self.set_child(y, k_left, z_left);
            self.set_parent(z_left, y);
            let z_color = self.get_color(z);
            self.set_color(y, z_color);
        }

        if y_was_black {
            self.delete_fixup(x, nil_parent);
        }
    }
}

/// HiGHS: HighsRbTree.h::highs::CacheMinRbTree<Impl> (lines 420-450)
///
/// Adds the cached minimum (`LinkType& first_`). C++ hides `RbTree::first()`,
/// `link(z, parent)` and `unlink(z)` with the cached versions, and the
/// inherited `link(z)` dispatches to the *derived* `link(z, parent)`.
///
/// The API enforces this. The base methods are named `rb_first` / `rb_link` /
/// `rb_link_parent` / `rb_unlink`, so the C++ names `first()`, `link(z)`,
/// `link_parent(z, parent)` and `unlink(z)` on a cache-min tree resolve only to
/// the cached versions below. The `cached_*` names are kept as aliases.
pub trait CacheMinRbTree: RbTree {
    /// C++: reading `LinkType& first_`
    fn first_cached_ref(&self) -> Self::LinkType;
    /// C++: writing `LinkType& first_`
    fn set_first_cached(&mut self, first: Self::LinkType);

    /// HiGHS: HighsRbTree.h::CacheMinRbTree::first (line 429)
    fn first(&self) -> Self::LinkType {
        self.first_cached_ref()
    }
    /// HiGHS: HighsRbTree.h::CacheMinRbTree::link (lines 432-442), `link(z, parent)`
    fn link_parent(&mut self, z: Self::LinkType, parent: Self::LinkType) {
        let no_link = Self::LinkType::NO_LINK;

        if self.first_cached_ref() == parent {
            if parent == no_link || PartialOrd::lt(&self.get_key(z), &self.get_key(parent)) {
                self.set_first_cached(z);
            }
        }

        self.rb_link_parent(z, parent);
    }
    /// HiGHS: HighsRbTree.h::RbTree::link(z) (lines 371-380) as seen from a CacheMinRbTree
    /// (`using RbTree<Impl>::link`; the final `Impl::link(z, y)` is the cached one).
    fn link(&mut self, z: Self::LinkType) {
        let no_link = Self::LinkType::NO_LINK;

        let mut y = no_link;
        let mut x = self.root_node();
        while x != no_link {
            y = x;
            let side = usize::from(PartialOrd::lt(&self.get_key(x), &self.get_key(z)));
            x = self.get_child(y, side);
        }

        CacheMinRbTree::link_parent(self, z, y);
    }
    /// HiGHS: HighsRbTree.h::CacheMinRbTree::unlink (lines 444-447)
    fn unlink(&mut self, z: Self::LinkType) {
        if z == self.first_cached_ref() {
            let first = self.first_cached_ref();
            let next = self.successor(first);
            self.set_first_cached(next);
        }
        self.rb_unlink(z);
    }

    /// Alias of [`CacheMinRbTree::first`].
    fn cached_first(&self) -> Self::LinkType {
        CacheMinRbTree::first(self)
    }
    /// Alias of [`CacheMinRbTree::link_parent`].
    fn cached_link_parent(&mut self, z: Self::LinkType, parent: Self::LinkType) {
        CacheMinRbTree::link_parent(self, z, parent);
    }
    /// Alias of [`CacheMinRbTree::link`].
    fn cached_link(&mut self, z: Self::LinkType) {
        CacheMinRbTree::link(self, z);
    }
    /// Alias of [`CacheMinRbTree::unlink`].
    fn cached_unlink(&mut self, z: Self::LinkType) {
        CacheMinRbTree::unlink(self, z);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test owner: node i has key `keys[i]`; links live in a vector.
    struct Tree {
        root: i32,
        first: i32,
        keys: Vec<i32>,
        links: Vec<RbTreeLinks<i32>>,
    }

    impl RbTree for Tree {
        type KeyType = i32;
        type LinkType = i32;
        fn root_node(&self) -> i32 {
            self.root
        }
        fn set_root_node(&mut self, root: i32) {
            self.root = root;
        }
        fn get_rb_tree_links(&self, node: i32) -> &RbTreeLinks<i32> {
            &self.links[node as usize]
        }
        fn get_rb_tree_links_mut(&mut self, node: i32) -> &mut RbTreeLinks<i32> {
            &mut self.links[node as usize]
        }
        fn get_key(&self, node: i32) -> i32 {
            self.keys[node as usize]
        }
    }

    impl CacheMinRbTree for Tree {
        fn first_cached_ref(&self) -> i32 {
            self.first
        }
        fn set_first_cached(&mut self, first: i32) {
            self.first = first;
        }
    }

    /// Checks parent links, no red-red edge and equal black height; returns the black height.
    fn check(t: &Tree, x: i32, parent: i32) -> usize {
        if x == -1 {
            return 1;
        }
        assert_eq!(t.get_parent(x), parent);
        let (l, r) = (t.get_child(x, 0), t.get_child(x, 1));
        if t.is_red(x) {
            assert!(t.is_black(l) && t.is_black(r));
        }
        let bl = check(t, l, x);
        let br = check(t, r, x);
        assert_eq!(bl, br);
        bl + usize::from(t.is_black(x))
    }

    #[test]
    fn links_parent_and_color_encoding() {
        let mut l = RbTreeLinks::<i32>::default();
        l.set_parent(-1);
        assert_eq!(l.get_parent(), -1);
        l.set_parent(41);
        assert_eq!(l.get_parent(), 41);
        l.make_red();
        assert!(l.is_red());
        assert_eq!(l.get_parent(), 41);
        l.set_color(false);
        assert!(l.is_black());
        assert_eq!(l.get_parent(), 41);

        let mut l64 = RbTreeLinks::<i64>::default();
        l64.set_parent(-1);
        assert_eq!(l64.get_parent(), -1);
        l64.set_parent(1 << 40);
        l64.make_red();
        assert_eq!(l64.get_parent(), 1 << 40);
    }

    #[test]
    fn cached_insert_delete_matches_sorted_order() {
        let n = 300usize;
        // distinct keys: 7919 is coprime to 1009, so i * 7919 mod 1009 is distinct for i < 1009
        let keys: Vec<i32> = (0..n as i32).map(|i| (i * 7919) % 1009).collect();
        let mut sorted: Vec<i32> = (0..n as i32).collect();
        sorted.sort_by_key(|&i| keys[i as usize]);

        let mut t = Tree {
            root: -1,
            first: -1,
            keys: keys.clone(),
            links: vec![RbTreeLinks::default(); n],
        };
        for i in 0..n as i32 {
            t.cached_link(i);
        }
        check(&t, t.root, -1);
        assert!(t.is_black(t.root));

        let walk = |t: &Tree| {
            let mut out = Vec::new();
            let mut x = t.cached_first();
            while x != -1 {
                out.push(x);
                x = t.successor(x);
            }
            out
        };
        assert_eq!(walk(&t), sorted);
        assert_eq!(t.last(), *sorted.last().unwrap());
        for i in 0..n as i32 {
            assert_eq!(t.find(&keys[i as usize]), (i, true));
        }

        for k in 0..100 {
            let z = t.cached_first();
            assert_eq!(z, sorted[k]);
            t.cached_unlink(z);
            check(&t, t.root, -1);
        }
        assert_eq!(walk(&t), sorted[100..].to_vec());
        assert_eq!(t.cached_first(), sorted[100]);
    }

    /// The C++ names (`link`, `unlink`, `first`) on a cache-min tree keep the
    /// cached minimum equal to the structural minimum through a random mix of
    /// inserts and deletes of arbitrary nodes. Keys are distinct, as in every
    /// HiGHS cache-min tree (each key ends with a node/column index): with
    /// equal keys C++ `CacheMinRbTree::link` does not move `first_` (strict `<`),
    /// so the cache can lag the structure there, in C++ as here.
    #[test]
    fn cpp_names_keep_cached_first_in_sync() {
        let n = 400usize;
        let mut st: u64 = 0x2545_f491_4f6c_dd1d;
        let mut nxt = || {
            st ^= st << 13;
            st ^= st >> 7;
            st ^= st << 17;
            st
        };
        let keys: Vec<i32> = (0..n as i32).map(|i| (i * 7919) % 1009).collect();
        let mut t = Tree {
            root: -1,
            first: -1,
            keys: keys.clone(),
            links: vec![RbTreeLinks::default(); n],
        };
        let mut linked = vec![false; n];
        for _ in 0..6000 {
            let i = (nxt() % n as u64) as usize;
            if linked[i] {
                t.unlink(i as i32);
            } else {
                t.link(i as i32);
            }
            linked[i] = !linked[i];
            assert_eq!(t.first(), t.rb_first());
            check(&t, t.root, -1);
        }
        let mut walk_keys = Vec::new();
        let mut x = t.first();
        while x != -1 {
            walk_keys.push(keys[x as usize]);
            x = t.successor(x);
        }
        let mut model: Vec<i32> = (0..n).filter(|&i| linked[i]).map(|i| keys[i]).collect();
        model.sort();
        assert_eq!(walk_keys, model);
    }
}
