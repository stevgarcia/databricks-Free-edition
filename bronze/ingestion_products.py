import dlt

#producst expectations
product_rules={
  "rule1": "product_id is not null",
  "rule2": "price>0"}
 

@dlt.table(  name="products_stage")
@dlt.expect_all_or_drop(product_rules)

def products_stage():
  return spark.readStream.table("catalogotest.source.products")

