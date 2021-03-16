use std::fmt;
use std::fmt::{Debug, Display, Formatter};
use std::path::PathBuf;
use toml::value::{Array, Datetime, Table};
use toml::Value;

#[derive(Debug)]
pub enum ExpectedType {
    String,
    Float,
    Int,
    UInt,
    Table,
    Array,
    Bool,
    DateTime,
}

impl Display for ExpectedType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

pub enum ParseErr {
    KeyNotFound(String),
    UnexpectedString(ExpectedType, String),
    UnexpectedFloat(ExpectedType, f64),
    UnexpectedInt(ExpectedType, i64),
    UnexpectedTable(ExpectedType, Table),
    UnexpectedArray(ExpectedType, Array),
    UnexpectedNone(ExpectedType),
    UnexpectedBool(ExpectedType, bool),
    UnexpectedDateTime(ExpectedType, Datetime),
    UnknownExperiment(String),
    UnknownIntegrator(String),
    FileOpen(String),
    FileParse(String),
    NoSeed,
    TooManyMarkedVerts(usize, usize),
    UnknownRgtpDistribDefn(String),
}

impl Debug for ParseErr {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            ParseErr::KeyNotFound(k) => {
                write!(f, "key not found: {}", k)
            }
            ParseErr::UnexpectedString(e, v) => {
                write!(f, "expected {}, found String {}", e, v)
            }
            ParseErr::UnexpectedFloat(e, v) => {
                write!(f, "expected {}, found Float {}", e, v)
            }
            ParseErr::UnexpectedInt(e, v) => {
                write!(f, "expected {}, found Int {}", e, v)
            }
            ParseErr::UnexpectedTable(e, v) => {
                write!(f, "expected {}, found Table {:?}", e, v)
            }
            ParseErr::UnexpectedArray(e, v) => {
                write!(f, "expected {}, found Array {:?}", e, v)
            }
            ParseErr::UnexpectedNone(e) => {
                write!(f, "expected {}, found None", e)
            }
            ParseErr::UnexpectedBool(e, v) => {
                write!(f, "expected {}, found Bool {}", e, v)
            }
            ParseErr::UnexpectedDateTime(e, v) => {
                write!(f, "expected {}, found DateTime {}", e, v)
            }
            ParseErr::UnknownExperiment(e) => {
                write!(f, "found unknown experiment: {}", e)
            }
            ParseErr::UnknownIntegrator(i) => {
                write!(f, "found unknown integrator: {}", i)
            }
            ParseErr::FileOpen(s) => {
                write!(f, "could not open file: {}", s)
            }
            ParseErr::FileParse(s) => {
                write!(f, "error parsing file as toml: {}", s)
            }
            ParseErr::NoSeed => {
                write!(
                    f,
                    "need a minimum of one seed, but no seeds found"
                )
            }
            ParseErr::TooManyMarkedVerts(max, found) => {
                write!(
                    f,
                    "too many marked verts. NVERTS={}, found={}",
                    max, found
                )
            }
            ParseErr::UnknownRgtpDistribDefn(s) => {
                write!(
                    f,
                    "unknown RhoGTPase distribution definition: {}",
                    s
                )
            }
        }
    }
}

impl ParseErr {
    fn unexpected(expected: ExpectedType, v: &Value) -> ParseErr {
        match v {
            Value::String(v) => {
                ParseErr::UnexpectedString(expected, v.clone())
            }
            Value::Integer(v) => {
                ParseErr::UnexpectedInt(expected, *v)
            }
            Value::Float(v) => {
                ParseErr::UnexpectedFloat(expected, *v)
            }
            Value::Boolean(v) => {
                ParseErr::UnexpectedBool(expected, *v)
            }
            Value::Datetime(v) => {
                ParseErr::UnexpectedDateTime(expected, v.clone())
            }
            Value::Array(v) => {
                ParseErr::UnexpectedArray(expected, v.clone())
            }
            Value::Table(v) => {
                ParseErr::UnexpectedTable(expected, v.clone())
            }
        }
    }
}

pub fn get_value(
    raw_table: &Table,
    key: &str,
) -> Result<Value, ParseErr> {
    if raw_table.contains_key(key) {
        if let Some(v) = raw_table.get(key) {
            Ok(v.clone())
        } else {
            Err(ParseErr::UnexpectedNone(ExpectedType::Table))
        }
    } else {
        Err(ParseErr::KeyNotFound(key.into()))
    }
}

pub fn parse_string(value: &Value) -> Result<String, ParseErr> {
    match value.as_str() {
        Some(s) => Ok(s.into()),
        None => {
            Err(ParseErr::unexpected(ExpectedType::String, value))
        }
    }
}

pub fn parse_f64(value: &Value) -> Result<f64, ParseErr> {
    match value.as_float() {
        Some(inner) => Ok(inner),
        None => parse_i64(value).map(|v| v as f64),
    }
}

pub fn parse_i64(value: &Value) -> Result<i64, ParseErr> {
    match value.as_integer() {
        Some(inner) => Ok(inner),
        None => Err(ParseErr::unexpected(ExpectedType::Int, value)),
    }
}

pub fn parse_u64(value: &Value) -> Result<u64, ParseErr> {
    match parse_i64(value) {
        Ok(inner) => {
            if inner < 0 {
                Err(ParseErr::UnexpectedInt(
                    ExpectedType::UInt,
                    inner,
                ))
            } else {
                Ok(inner as u64)
            }
        }
        Err(e) => Err(e),
    }
}

pub fn parse_bool(value: &Value) -> Result<bool, ParseErr> {
    match value.as_bool() {
        Some(inner) => Ok(inner),
        None => Err(ParseErr::unexpected(ExpectedType::Bool, value)),
    }
}

pub fn parse_table(value: &Value) -> Result<Table, ParseErr> {
    match value.as_table() {
        Some(inner) => Ok(inner.clone()),
        None => Err(ParseErr::unexpected(ExpectedType::Table, value)),
    }
}

pub fn parse_array(value: &Value) -> Result<Vec<Value>, ParseErr> {
    match value.as_array() {
        Some(inner) => Ok(inner.clone()),
        None => Err(ParseErr::unexpected(ExpectedType::Array, value)),
    }
}

pub fn parse_path(value: &Value) -> Result<PathBuf, ParseErr> {
    parse_string(value).map(|s| PathBuf::from(s))
}

pub fn get_optional_f64(
    raw_table: &Table,
    key: &str,
) -> Result<Option<f64>, ParseErr> {
    match get_value(raw_table, key) {
        Err(ParseErr::KeyNotFound(_)) => Ok(None),
        Err(e) => Err(e),
        Ok(v) => parse_f64(&v).map(Some),
    }
}

pub fn get_path(
    raw_table: &Table,
    key: &str,
) -> Result<Option<PathBuf>, ParseErr> {
    match get_value(raw_table, key) {
        Err(ParseErr::KeyNotFound(_)) => Ok(None),
        Err(e) => Err(e),
        Ok(v) => parse_string(&v).map(|s| Some(PathBuf::from(s))),
    }
}

pub fn get_optional_path(
    raw_table: &Table,
    key: &str,
) -> Result<Option<PathBuf>, ParseErr> {
    match get_value(raw_table, key) {
        Err(ParseErr::KeyNotFound(_)) => Ok(None),
        Err(e) => Err(e),
        Ok(v) => {
            let p = parse_path(&v)?;
            Ok(Some(p))
        }
    }
}

pub fn get_optional_usize(
    raw_table: &Table,
    key: &str,
) -> Result<Option<usize>, ParseErr> {
    match get_value(raw_table, key) {
        Err(ParseErr::KeyNotFound(_)) => Ok(None),
        Err(e) => Err(e),
        Ok(v) => {
            let i = parse_i64(&v)?;
            if i < 0 {
                Err(ParseErr::UnexpectedInt(ExpectedType::UInt, i))
            } else {
                Ok(Some(i as usize))
            }
        }
    }
}

pub fn get_optional_bool(
    raw_table: &Table,
    key: &str,
) -> Result<Option<bool>, ParseErr> {
    match get_value(raw_table, key) {
        Err(ParseErr::KeyNotFound(_)) => Ok(None),
        Err(e) => Err(e),
        Ok(v) => {
            let b = parse_bool(&v)?;
            Ok(Some(b))
        }
    }
}
