use crate::utils::pcg32::Pcg32;
use rand::distributions::Distribution;
use rand_distr::Normal;
use serde::de::Visitor;
use serde::{de, Deserialize, Deserializer, Serialize};
use std::fmt;

#[derive(Clone, Copy, Serialize)]
pub struct NormalDistrib {
    mean: f64,
    std: f64,
    #[serde(skip)]
    distrib: Normal<f64>,
}

impl NormalDistrib {
    pub fn sample(&self, rng: &mut Pcg32) -> f64 {
        self.distrib.sample(rng)
    }

    pub fn new(mean: f64, std: f64) -> NormalDistrib {
        NormalDistrib {
            mean,
            std,
            distrib: Normal::new(mean, std).unwrap(),
        }
    }
}

impl<'de> Deserialize<'de> for NormalDistrib {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        #[serde(field_identifier, rename_all = "lowercase")]
        enum Field {
            Mean,
            Std,
        }

        struct NdVisitor;

        impl<'de> Visitor<'de> for NdVisitor {
            type Value = NormalDistrib;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("struct NormalDistrib")
            }

            fn visit_seq<V>(self, mut seq: V) -> Result<NormalDistrib, V::Error>
            where
                V: de::SeqAccess<'de>,
            {
                let mean = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let std = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(1, &self))?;
                Ok(NormalDistrib {
                    mean,
                    std,
                    distrib: Normal::new(mean, std).unwrap(),
                })
            }

            fn visit_map<V>(self, mut map: V) -> Result<NormalDistrib, V::Error>
            where
                V: de::MapAccess<'de>,
            {
                let mut mean = None;
                let mut std = None;
                while let Some(key) = map.next_key()? {
                    match key {
                        Field::Mean => {
                            if mean.is_some() {
                                return Err(de::Error::duplicate_field("mean"));
                            }
                            mean = Some(map.next_value()?);
                        }
                        Field::Std => {
                            if std.is_some() {
                                return Err(de::Error::duplicate_field("std"));
                            }
                            std = Some(map.next_value()?);
                        }
                    }
                }
                let mean =
                    mean.ok_or_else(|| de::Error::missing_field("mean"))?;
                let std = std.ok_or_else(|| de::Error::missing_field("std"))?;
                Ok(NormalDistrib {
                    mean,
                    std,
                    distrib: Normal::new(mean, std).unwrap(),
                })
            }
        }

        const FIELDS: &[&str] = &["mean", "std"];
        deserializer.deserialize_struct("NormalDistrib", FIELDS, NdVisitor)
    }
}

impl PartialEq for NormalDistrib {
    fn eq(&self, other: &Self) -> bool {
        self.mean == other.mean && self.std == other.std
    }
}
