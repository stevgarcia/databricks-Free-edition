import dlt
from  pyspark.sql.functions import *


@dlt.table(name="business_sales")

def business_sales():
  
    df_fact=spark.read.table("Fact_sales")
    df_dim=spark.read.table("dim_customers")
    df_dimProd=spark.read.table("dim_products")

    df_join=df_fact.join(df_dim, df_fact.customer_id == df_dim.customer_id).join(df_dimProd, df_fact.product_id == df_dimProd.product_id)

    df_prun=df_join.select("region", "category", "total_amount")

    

    df_agg=df_prun.groupBy("region", "category").agg(sum("total_amount").alias("total_sales"))

    return df_agg