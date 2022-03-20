import os, logging
from datetime import date

def dump_database(MONGO_INITDB_ROOT_USERNAME, MONGO_INITDB_ROOT_PASSWORD):
	"""
	function for dumping the database. needs root access
	"""
	today = date.today()
	# dd.mm.YY
	current_date = today.strftime("%d.%m.%Y")
	path = os.getcwd()
	path = os.path.join(path, "mongodumps", "{}".format(current_date))
	logging.info(path)
	# exec_into_docker_command = "docker exec -it mongodb bash"
	# logging.info(os.system(exec_into_docker_command))

	# запуск команды по дампу
	command = "mongodump --host {} -u {} -p {} --authenticationDatabase admin -o={}".format("mongodb_api", MONGO_INITDB_ROOT_USERNAME, MONGO_INITDB_ROOT_PASSWORD, path)
	response = os.system(command)
	logging.info(response)

	# запуск команды архивирования
	command = f"zip -r {path}.zip {path}"
	response = os.system(command)
	logging.info(response)

	return (response, path)