# Databricks delta live tables and pipelines

Welcome to the ##Databricks project##

This project demonstrates a comprehensive ETL solution, it leverages its data with delta live tables, expectations, AUTOCDC and incremental loads.


<img width="955" height="592" alt="image" src="https://github.com/user-attachments/assets/22e62915-698d-41cb-bd34-eb01e4efa20f" />

## 📖 Project Overview
This project involves:

1. **Unit catalog**: a Data warehouse created in Unity catalog over the source schema.
2. **ETL Pipeline**: a pipeline defined in databricks with multiple notebooks for creating a medallion architecture.
3. **Data warehouse**: a predefined SQL script that creates a east sales, west sales, customers and products


   
#### Objective
Build a databricks pipeline to support an ETL with a medallion architecture and the use of Delta live tables .

#### Specifications
- **Data Sources**: a predefined SQL script that creates a east sales, west sales, customers and products in UNITY CATALOG.

<img width="726" height="628" alt="image" src="https://github.com/user-attachments/assets/bb029c14-9d9d-4ea3-92f7-81f37d13bfcf" />

- **Bronze layer**: involves the use of an empty streaming table and an append flow to obtain the east sales and west sales, and unify them. 
- It creates a streaming table to get the products data with defined EXPECTATIONS.
- It creates a streaming table to get the customers data with defined EXPECTATIONS.

- **Silver layer**: it creates an enriched streamings views on top of the streaming tables from the bronze layer. It supports the Change data capture with the use of the flow AUTOCDC and apply simple transformations.
