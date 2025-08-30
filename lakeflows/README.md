# Databricks lakeFlow Jobs

Welcome to the ##lakeFlow Jobs##

This project creates a dataflow job, where two task are connected, the first task Selects with SQL the name of multiple files and with the use of parameters the second one, gets the names with the use of widgets and then using Pyspark writes the data in delta format in a volume as a sink.


<img width="823" height="275" alt="image" src="https://github.com/user-attachments/assets/033e04ad-7ec0-44f6-84e9-6a04cf8c9e59" />



## 📖 Project Overview
This project involves:

1. **Volume**: his project has 3 parquet files, orders, products and regions hosted in a volume named as  source.
2.  **Notebook**: a Pyspark notebook named ingestion that uses widgets expecting the parameter file_names and then writes all the parquet files over a volume in a folder called sink.

<img width="644" height="470" alt="image" src="https://github.com/user-attachments/assets/336b4b1d-fe65-44ee-9675-284d35017c29" />


3. **SQL mapping**: an SQL scrip that creates a table in the default schema with the file names
   
<img width="1028" height="712" alt="image" src="https://github.com/user-attachments/assets/a8078202-d555-4002-98ec-6b8ea8a364cd" />

4. **SQL task**: the task executes an select * from mapping
5. **Loop**: a loop  named DynamicIngestion that iterates with {{tasks.SQLArray.output.rows}}
6.  **task inside a loop**: a task named iteration that sends the file names to a Notebook previously configured with widgets to load the data in the sink folder in delta format. 


<img width="1037" height="516" alt="image" src="https://github.com/user-attachments/assets/16ea6293-f7ee-46e6-ad84-b49afa7583cb" />





