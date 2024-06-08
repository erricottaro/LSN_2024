import subprocess
import numpy as np

#avvia la simulazione a temperatura T
def avvia_simulatore():
    # Comando Bash per eseguire il programma es08.1.x
    comando_bash = "./es08.1.x"
    
    # Esegue il comando Bash
    print("Running simulation...")
    subprocess.run(comando_bash, shell=True)
    print("Done")

#Copia i file di output in apposite cartelle
def copia_output():
    dir_name = f"OPTIMIZED"
    comando_dir = "mkdir -p "+dir_name
    comando_copia = "cp OUTPUT/distribution.dat OUTPUT/energy.dat OUTPUT/acceptance.dat OUTPUT/output.dat "+dir_name
    print("Building directory...")
    subprocess.run(comando_dir, shell=True)
    print("Done")
    print("Copying output...")
    subprocess.run(comando_copia, shell=True)
    print("Done")


if __name__ == "__main__":

    avvia_simulatore()
    copia_output()
