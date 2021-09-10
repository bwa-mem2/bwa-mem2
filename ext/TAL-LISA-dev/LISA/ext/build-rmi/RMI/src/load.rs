// < begin copyright > 
// Copyright Ryan Marcus 2020
// 
// See root directory of this project for license terms.
// 
// < end copyright > 
 
 
use memmap::MmapOptions;
use rmi_lib::{RMITrainingData, RMITrainingDataIteratorProvider, KeyType};
use byteorder::{LittleEndian, ReadBytesExt};
use std::fs::File;
use std::convert::TryInto;

pub enum DataType {
    UINT64,
    UINT32,
    FLOAT64
}

struct SliceAdapterU64 {
    data: memmap::Mmap,
    length: usize
}

impl RMITrainingDataIteratorProvider for SliceAdapterU64 {
    type InpType = u64;
    fn cdf_iter(&self) -> Box<dyn Iterator<Item = (Self::InpType, usize)> + '_> {
        Box::new((0..self.length).map(move |i| self.get(i).unwrap()))
    }
    
    fn get(&self, idx: usize) -> Option<(Self::InpType, usize)> {
        if idx >= self.length { return None; };
        let mi = u64::from_le_bytes((&self.data[8 + idx * 8..8 + (idx + 1) * 8])
                                    .try_into().unwrap());
        return Some((mi.into(), idx));
    }
    
    fn key_type(&self) -> KeyType {
        KeyType::U64
    }
    
    fn len(&self) -> usize { self.length }
}


struct SliceAdapterU32 {
    data: memmap::Mmap,
    length: usize
}

impl RMITrainingDataIteratorProvider for SliceAdapterU32 {
    type InpType = u32;
    fn cdf_iter(&self) -> Box<dyn Iterator<Item = (Self::InpType, usize)> + '_> {
        Box::new((0..self.length).map(move |i| self.get(i).unwrap()))
    }
    
    fn get(&self, idx: usize) -> Option<(Self::InpType, usize)> {
        if idx >= self.length { return None; };
        let mi = (&self.data[8 + idx * 4..8 + (idx + 1) * 4])
            .read_u32::<LittleEndian>().unwrap().into();
        return Some((mi, idx));
    }
    
    fn key_type(&self) -> KeyType {
        KeyType::U32
    }
    
    fn len(&self) -> usize { self.length }
}

struct SliceAdapterF64 {
    data: memmap::Mmap,
    length: usize
}

impl RMITrainingDataIteratorProvider for SliceAdapterF64 {
    type InpType = f64;
    fn cdf_iter(&self) -> Box<dyn Iterator<Item = (Self::InpType, usize)> + '_> {
        Box::new((0..self.length).map(move |i| self.get(i).unwrap()))
    }
    
    fn get(&self, idx: usize) -> Option<(Self::InpType, usize)> {
        if idx >= self.length { return None; };
        let mi = (&self.data[8 + idx * 8..8 + (idx + 1) * 8])
            .read_f64::<LittleEndian>().unwrap().into();
        return Some((mi, idx));
    }
    
    fn key_type(&self) -> KeyType {
        KeyType::F64
    }
    
    fn len(&self) -> usize { self.length }
}

pub enum RMIMMap {
    UINT64(RMITrainingData<u64>),
    UINT32(RMITrainingData<u32>),
    FLOAT64(RMITrainingData<f64>)
}

macro_rules! dynamic {
    ($funcname: expr, $data: expr $(, $p: expr )*) => {
        match $data {
            load::RMIMMap::UINT64(mut x) => $funcname(&mut x, $($p),*),
            load::RMIMMap::UINT32(mut x) => $funcname(&mut x, $($p),*),
            load::RMIMMap::FLOAT64(mut x) => $funcname(&mut x, $($p),*),
        }
    }
}


impl RMIMMap {
    pub fn soft_copy(&self) -> RMIMMap {
        match self {
            RMIMMap::UINT64(x) => RMIMMap::UINT64(x.soft_copy()),
            RMIMMap::UINT32(x) => RMIMMap::UINT32(x.soft_copy()),
            RMIMMap::FLOAT64(x) => RMIMMap::FLOAT64(x.soft_copy()),
        }
    }

    pub fn into_u64(self) -> Option<RMITrainingData<u64>> {
        match self {
            RMIMMap::UINT64(x) => Some(x),
            _ => None
        }
    }
}
                

pub fn load_data(filepath: &str,
                 dt: DataType) -> (usize, RMIMMap) {
    let fd = File::open(filepath).unwrap_or_else(|_| {
        panic!("Unable to open data file at {}", filepath)
    });

    let mmap = unsafe { MmapOptions::new().map(&fd).unwrap() };
    let num_items = (&mmap[0..8]).read_u64::<LittleEndian>().unwrap() as usize;

    let rtd = match dt {
        DataType::UINT64 =>
            RMIMMap::UINT64(RMITrainingData::new(Box::new(
                SliceAdapterU64 { data: mmap, length: num_items }
            ))),
        DataType::UINT32 =>
            RMIMMap::UINT32(RMITrainingData::new(Box::new(
                SliceAdapterU32 { data: mmap, length: num_items }
            ))),
        DataType::FLOAT64 =>
            RMIMMap::FLOAT64(RMITrainingData::new(Box::new(
                SliceAdapterF64 { data: mmap, length: num_items }
            )))
    };

    return (num_items, rtd);
}
