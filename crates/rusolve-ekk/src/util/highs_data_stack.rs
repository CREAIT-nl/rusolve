//! HiGHS source: util/HighsDataStack.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The byte stack that `HighsPostsolveStack` serialises its reductions into.
//! Complete in the skeleton (no fill-in needed).
//!
//! C++ `push(const T&)` `memcpy`s a trivially copyable `T`; `pop(T&)` moves
//! `position` back by `sizeof(T)` and copies the bytes out. The Rust port
//! keeps the *same stack discipline* through the [`HighsDataStackItem`]
//! trait: every item has a fixed byte size `SIZE`, `push` appends its bytes,
//! `pop` moves `position` back by `SIZE` and reads forward from there.
//! Structs are serialised field by field in declaration order (no padding),
//! so byte offsets differ from the C++ `sizeof` layout. This is unobservable:
//! positions are only ever produced by `get_current_data_size` and consumed by
//! `set_position`/`pop` of the same item sequence (unfaithful by design).
//!
//! `push(const std::vector<T>&)` / `pop(std::vector<T>&)` are
//! [`HighsDataStack::push_vec`] / [`HighsDataStack::pop_vec`]: the elements,
//! then the element count as a `usize`; popping a count of 0 clears the vector.

use crate::prelude::*;

/// A value that can be pushed onto / popped from a [`HighsDataStack`]
/// (C++: any `IS_TRIVIALLY_COPYABLE(T)`).
pub trait HighsDataStackItem: Sized + Default + Clone {
    /// Number of bytes the item occupies on the stack.
    const SIZE: usize;
    /// Append the bytes of `self` to `data`.
    fn write_bytes(&self, data: &mut Vec<u8>);
    /// Read one item forward from `data[*p..]`, advancing `*p` by `SIZE`.
    fn read_bytes(data: &[u8], p: &mut usize) -> Self;
}

macro_rules! impl_data_stack_item_num {
    ($($t:ty),*) => {
        $(
            impl HighsDataStackItem for $t {
                const SIZE: usize = std::mem::size_of::<$t>();
                fn write_bytes(&self, data: &mut Vec<u8>) {
                    data.extend_from_slice(&self.to_ne_bytes());
                }
                fn read_bytes(data: &[u8], p: &mut usize) -> Self {
                    let mut bytes = [0u8; std::mem::size_of::<$t>()];
                    bytes.copy_from_slice(&data[*p..*p + Self::SIZE]);
                    *p += Self::SIZE;
                    <$t>::from_ne_bytes(bytes)
                }
            }
        )*
    };
}
impl_data_stack_item_num!(f64, i32, u8, u32, u64, usize, i64);

impl HighsDataStackItem for bool {
    const SIZE: usize = 1;
    fn write_bytes(&self, data: &mut Vec<u8>) {
        data.push(u8::from(*self));
    }
    fn read_bytes(data: &[u8], p: &mut usize) -> Self {
        let v = data[*p] != 0;
        *p += 1;
        v
    }
}

impl HighsDataStackItem for HighsBasisStatus {
    const SIZE: usize = 1;
    fn write_bytes(&self, data: &mut Vec<u8>) {
        data.push(*self as u8);
    }
    fn read_bytes(data: &[u8], p: &mut usize) -> Self {
        let v = match data[*p] {
            0 => HighsBasisStatus::kLower,
            1 => HighsBasisStatus::kBasic,
            2 => HighsBasisStatus::kUpper,
            3 => HighsBasisStatus::kZero,
            _ => HighsBasisStatus::kNonbasic,
        };
        *p += 1;
        v
    }
}

/// Implements [`HighsDataStackItem`] for a plain struct, serialising the
/// listed fields in the given (declaration) order:
/// `impl_highs_data_stack_item!(Type { field: FieldType, ... });`
macro_rules! impl_highs_data_stack_item {
    ($t:ty { $($f:ident : $ft:ty),* $(,)? }) => {
        impl $crate::util::highs_data_stack::HighsDataStackItem for $t {
            const SIZE: usize =
                0 $(+ <$ft as $crate::util::highs_data_stack::HighsDataStackItem>::SIZE)*;
            fn write_bytes(&self, data: &mut Vec<u8>) {
                $(
                    $crate::util::highs_data_stack::HighsDataStackItem::write_bytes(
                        &self.$f, data,
                    );
                )*
            }
            fn read_bytes(data: &[u8], p: &mut usize) -> Self {
                // Struct literal fields are evaluated in source order.
                Self {
                    $(
                        $f: <$ft as $crate::util::highs_data_stack::HighsDataStackItem>::read_bytes(
                            data, p,
                        ),
                    )*
                }
            }
        }
    };
}
pub(crate) use impl_highs_data_stack_item;

/// HiGHS: HighsDataStack.h::HighsDataStack
///
/// C++ leaves `position` uninitialised; it is always set by `resetPosition` /
/// `setPosition` before a `pop`, so the Rust default `0` is unobservable.
#[derive(Clone, Debug, Default)]
pub struct HighsDataStack {
    data: Vec<u8>,
    position: usize,
}

impl HighsDataStack {
    /// HiGHS: HighsDataStack.h::HighsDataStack::resetPosition
    pub fn reset_position(&mut self) {
        self.position = self.data.len();
    }

    /// HiGHS: HighsDataStack.h::HighsDataStack::push(const T&)
    pub fn push<T: HighsDataStackItem>(&mut self, r: &T) {
        r.write_bytes(&mut self.data);
    }

    /// HiGHS: HighsDataStack.h::HighsDataStack::pop(T&)
    pub fn pop<T: HighsDataStackItem>(&mut self, r: &mut T) {
        self.position -= T::SIZE;
        let mut p = self.position;
        *r = T::read_bytes(&self.data, &mut p);
    }

    /// HiGHS: HighsDataStack.h::HighsDataStack::push(const std::vector<T>&)
    pub fn push_vec<T: HighsDataStackItem>(&mut self, r: &[T]) {
        for x in r {
            x.write_bytes(&mut self.data);
        }
        r.len().write_bytes(&mut self.data);
    }

    /// HiGHS: HighsDataStack.h::HighsDataStack::pop(std::vector<T>&)
    pub fn pop_vec<T: HighsDataStackItem>(&mut self, r: &mut Vec<T>) {
        self.position -= usize::SIZE;
        let mut p = self.position;
        let num_data = usize::read_bytes(&self.data, &mut p);
        if num_data == 0 {
            r.clear();
        } else {
            r.clear();
            self.position -= num_data * T::SIZE;
            let mut p = self.position;
            r.reserve(num_data);
            for _ in 0..num_data {
                r.push(T::read_bytes(&self.data, &mut p));
            }
        }
    }

    /// HiGHS: HighsDataStack.h::HighsDataStack::setPosition
    pub fn set_position(&mut self, position_: usize) {
        self.position = position_;
    }

    /// HiGHS: HighsDataStack.h::HighsDataStack::getCurrentDataSize
    pub fn get_current_data_size(&self) -> usize {
        self.data.len()
    }
}
