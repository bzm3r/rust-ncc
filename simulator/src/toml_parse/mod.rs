use std::path::PathBuf;
use toml::value::{Array, Datetime, Table};
use toml::Value;

#[derive(Debug)]
pub enum ParseErr {
    KeyNotFound(String),
    UnexpectedString(String),
    UnexpectedFloat(f64),
    UnexpectedInt(i64),
    UnexpectedTable(Table),
    UnexpectedArray(Array),
    UnexpectedNone,
    UnexpectedBool(bool),
    UnexpectedDateTime(Datetime),
    UnknownExperiment(String),
    UnknownIntegrator(String),
    FileOpen(String),
    FileParse(String),
    MultipleNoneSeeds,
}

impl ParseErr {
    fn unexpected(v: &Value) -> ParseErr {
        match v {
            Value::String(v) => ParseErr::UnexpectedString(v.clone()),
            Value::Integer(v) => ParseErr::UnexpectedInt(*v),
            Value::Float(v) => ParseErr::UnexpectedFloat(*v),
            Value::Boolean(v) => ParseErr::UnexpectedBool(*v),
            Value::Datetime(v) => {
                ParseErr::UnexpectedDateTime(v.clone())
            }
            Value::Array(v) => ParseErr::UnexpectedArray(v.clone()),
            Value::Table(v) => ParseErr::UnexpectedTable(v.clone()),
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
            Err(ParseErr::UnexpectedNone)
        }
    } else {
        Err(ParseErr::KeyNotFound(key.into()))
    }
}

pub fn parse_string(value: &Value) -> Result<String, ParseErr> {
    match value.as_str() {
        Some(s) => Ok(s.into()),
        None => Err(ParseErr::unexpected(value)),
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
        None => Err(ParseErr::unexpected(value)),
    }
}

pub fn parse_u64(value: &Value) -> Result<u64, ParseErr> {
    match parse_i64(value) {
        Ok(inner) => {
            if inner < 0 {
                Err(ParseErr::UnexpectedInt(inner))
            } else {
                Ok(inner as u64)
            }
        }
        Err(e) => Err(e),
    }
}

pub fn parse_table(value: &Value) -> Result<Table, ParseErr> {
    match value.as_table() {
        Some(inner) => Ok(inner.clone()),
        None => Err(ParseErr::unexpected(value)),
    }
}

pub fn parse_array(value: &Value) -> Result<Vec<Value>, ParseErr> {
    match value.as_array() {
        Some(inner) => Ok(inner.clone()),
        None => Err(ParseErr::unexpected(value)),
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
                Err(ParseErr::UnexpectedInt(i))
            } else {
                Ok(Some(i as usize))
            }
        }
    }
}
