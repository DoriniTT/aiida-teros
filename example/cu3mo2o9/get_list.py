from aiida.orm import List, load_node
from aiida import load_profile

# Load the AiiDA profile (necessary to interact with the AiiDA database)
load_profile()

# Create an AiiDA List
aiida_list = load_node(223760)

# Store the AiiDA List in the database
aiida_list.store()

# Retrieve the list as a Python list
python_list = aiida_list.get_list()

print(python_list)