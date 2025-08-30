# Databricks lakeFlow Jobs

Welcome to the ##lakeFlow Jobs##

This project has 3 parquet files, orders, products and regions hosted in a volume called source.





## 📖 Project Overview
This project involves:

1. **Volume**: his project has 3 parquet files, orders, products and regions hosted in a volume named as  source.
2.  **Notebook**: a Pyspark notebook named ingestion that uses widgets expecting the parameter file_names and then writes all the parquet files over a volume in a folder called sink.

<img width="644" height="470" alt="image" src="https://github.com/user-attachments/assets/336b4b1d-fe65-44ee-9675-284d35017c29" />


3. **SQL mapping**: an SQL scrip that creates a table in the default schema with the file names
   
<img width="1028" height="712" alt="image" src="https://github.com/user-attachments/assets/a8078202-d555-4002-98ec-6b8ea8a364cd" />

4. **SQL task**: the task executes an select * from mapping
5. **Loop**: a loop  named DynamicIngestion that iterates with {{tasks.SQLArray.output.rows}}
6.  **task inside a loop**: a task named iteration that sends the file names to a Notebook previously configured with widgets to load the data in the sink folder. 

   
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

