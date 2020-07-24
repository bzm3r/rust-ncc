import fastavro
import json

test_recs = []
with open('test_schema.avsc') as sf:
    fastavro.parse_schema(json.load(sf))
    with open('test_dat.avro', 'rb') as df:
        for r in fastavro.reader(df):
            test_recs.append(r)
