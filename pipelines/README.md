# Databricks delta live tables and pipelines

Welcome to the ##Databricks project##

This project demonstrates a comprehensive ETL solution, it leverages its data with delta live tables, allowing incremental loads.


<img width="955" height="592" alt="image" src="https://github.com/user-attachments/assets/22e62915-698d-41cb-bd34-eb01e4efa20f" />

## 📖 Project Overview
This project involves:

1. **Unit catalog**: Defining a catalog with a schema named source.
2. **ETL Pipeline**: a pipeline defined in databricks with multiple notebooks for creating a medallion architecture.
3. **Data warehouse**: a predefined SQL script that creates a east sales, west sales, customers and products
4. **Bronze layer**: where the sales from east and west are unified, the product data and customers data.
5. **Silver layer**: where we apply some basic transformations using pyspark and upserts.
6. **Gold layer**: where we apply some basic transformations using pyspark and upserts

   
#### Objective
Build a databricks pipeline to .

#### Specifications
- **Data Sources**: a predefined SQL script that creates a east sales, west sales, customers and products in UNITY CATALOG.

<img width="726" height="628" alt="image" src="https://github.com/user-attachments/assets/bb029c14-9d9d-4ea3-92f7-81f37d13bfcf" />

- **Data Quality**: Cleanse and resolve data quality issues prior to analysis.
- **Integration**: Combine both sources into a single, user-friendly data model designed for analytical queries.
- **Scope**: Focus on the latest dataset only; historization of data is not required.
- **Documentation**: Provide clear documentation of the data model to support both business stakeholders and analytics teams.
