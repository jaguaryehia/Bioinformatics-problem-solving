import os
# path to the file to run in cmd line
cmd = "streamlit run --server.port 8080 C:\\Users\\jaguar\\PycharmProjects\\learningpandas\\AGraduationProject.py"

# returns the exit code in Windows
returned_value = os.system(cmd)
print('returned value:', returned_value)

