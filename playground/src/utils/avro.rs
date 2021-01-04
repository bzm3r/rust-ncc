use avro_rs::Schema;

pub trait Schematized {
    fn schematize(parent_space: Option<String>) -> Schema;
}

impl Schematized for bool {
    fn schematize(_parent_space: Option<String>) -> Schema {
        avro_rs::schema::Schema::Boolean
    }
}
