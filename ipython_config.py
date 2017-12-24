c = get_config()

#set matplotlib backend to nbagg
c.InteractiveShellApp.matplotlib = "nbagg"

c.InteractiveShellApp.exec_lines = [
"%pylab\n",
"%nbagg\n",
"import sygma as s\n",
"import omega as o\n",
"import nugridpy.nugridse as mp\n",
"import nugridpy.mesa as ms\n",
"from NuGrid_Mesa_Explorer_py3 import start_explorer\n",
"from SYGMA import start_SYGMA\n",
"from OMEGA import start_OMEGA\n",
"ms.set_nugrid_path('/data/nugrid_apod2')\n",
"mp.set_nugrid_path('/data/nugrid_apod2')\n",
"from IPython.display import HTML\n",
"display(HTML(\'<style>.container { width:80% !important; }</style>\'))\n"
]

c.NotebookNotary.data_dir="/home/user/.secret"
c.NotebookNotary.db_file= "/home/user/.secret/nbsignatures.db"
c.NotebookNotary.db_file= "/home/user/.secret/nbsignatures.db"
c.NotebookApp.secret_file="/home/user/.secret/notebook_secret"

#c.NotebookApp.secret_file="/home/user/.local/share/jupyter/notebook_secret"
#c.update({"load_extensions":{"init_cell/main":True}})
