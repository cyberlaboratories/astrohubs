c = get_config()

#set matplotlib backend to nbagg          
c.InteractiveShellApp.matplotlib = "nbagg"

c.InteractiveShellApp.exec_lines = [
"%pylab inline\n",
"%nbagg\n",
"import sygma as s\n",
"import omega as o\n",
"from SYGMA import start_SYGMA\n",
]

c.NotebookNotary.data_dir="/home/user/.secret"
c.NotebookNotary.db_file= "/home/user/.secret/nbsignatures.db"
c.NotebookNotary.db_file= "/home/user/.secret/nbsignatures.db"
c.NotebookApp.secret_file="/home/user/.secret/notebook_secret"

#c.NotebookApp.secret_file="/home/user/.local/share/jupyter/notebook_secret"    
#c.update({"load_extensions":{"init_cell/main":True}}) 
