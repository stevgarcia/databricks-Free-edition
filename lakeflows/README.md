# Databricks lakeFlow Jobs

Welcome to the ##lakeFlow Jobs##

This project demonstrates 




## 📖 Project Overview
This project involves:

1. **Unit catalog**: a Data warehouse created in Unity catalog over the source schema.
2. **ETL Pipeline**: a pipeline defined in databricks with multiple notebooks for creating a medallion architecture.
3. **Data warehouse**: a predefined SQL script that creates a east sales, west sales, customers and products


   
#### Objective
Build a databricks pipeline to support an ETL with a medallion architecture and the use of Delta live tables .




#### Specifications
- **Data Sources**: a predefined SQL script that creates a east sales, west sales, customers and products in UNITY CATALOG.



- **Bronze layer**: involves the use of an empty streaming table and an append flow to obtain the east sales and west sales, and unify them. 
- It creates a streaming table to get the products data with defined EXPECTATIONS.
- It creates a streaming table to get the customers data with defined EXPECTATIONS.

- **Silver layer**: it creates an enriched streamings views on top of the streaming tables from the bronze layer. It supports the Change data capture with the use of the flow AUTOCDC and applies simple transformations.

- **Gold layer**: it creates the dimensions with the user of streaming tables on top of the enriched views from silver layer usin AUTOCDC and supporting SCD type 2.

- **Business view**: the business_sales view is a materialized view that comes from the join etween fact_sales and dim_product

