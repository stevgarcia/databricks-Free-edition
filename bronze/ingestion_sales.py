import dlt

sales_rules= {
    "rule1":"sales_id is not null"
    
}

#empty streaming table
dlt.create_streaming_table(  name="append_sales", expect_all_or_drop=sales_rules)



#creating east sales flow

@dlt.append_flow(target="append_sales")
def east_sales():
  return spark.readStream.table("catalogotest.source.sales_east")

#creating awest sales flow

@dlt.append_flow(target="append_sales")
def west_sales():
  return spark.readStream.table("catalogotest.source.sales_east")