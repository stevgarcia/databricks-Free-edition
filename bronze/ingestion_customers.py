import dlt

#customers expectations
customer_rules={
  "rule1": "customer_id is not null",
  "rule2": "customer_name is not null"}
 

@dlt.table(  name="customers_stage")
@dlt.expect_all(customer_rules)
def customers_stage():
  return spark.readStream.table("catalogotest.source.customers")