import subprocess
import numpy as np

#avvia la simulazione a temperatura T
def avvia_simulatore():
    # Comando Bash per eseguire il programma simulator.exe
    comando_bash = "./simulator.exe"
    
    # Esegue il comando Bash
    print("Running simulation...")
    subprocess.run(comando_bash, shell=True)
    print("Done")

#Modifica la temperatura nel file di input
def modifica_valore_input(filename, nuovo_valore):
    with open(filename, 'r') as file:
        linee = file.readlines()

    with open(filename, 'w') as file:
        for linea in linee:
            if linea.startswith("TEMP"):
                # Modifica il valore numerico associato all'intestazione "T"
                nuova_linea = f"TEMP                   {nuovo_valore}\n"
                file.write(nuova_linea)
            else:
                # Scrivi le altre righe del file senza modificarle
                file.write(linea)

#Modifica l'opzione restart nel file di input
def modifica_restart(filename, value):
    with open(filename, 'r') as file:
        linee = file.readlines()

    with open(filename, 'w') as file:
        for linea in linee:
            if linea.startswith("RESTART"):
                # Modifica il valore numerico associato all'intestazione "RESTART"
                nuova_linea = f"RESTART                {value}\n"
                file.write(nuova_linea)
            else:
                # Scrivi le altre righe del file senza modificarle
                file.write(linea)

#Copia i file di output in apposite cartelle
def copia_output(T, method):
    dir_name = f"../../{method}"+f"{format(T, '.3f')}"
    comando_dir = "mkdir -p "+dir_name
    comando_copia = "cp ../OUTPUT/*.dat "+dir_name
    #print(comando_dir)
    #print(comando_copia)
    print("Building directory...")
    subprocess.run(comando_dir, shell=True)
    print("Done")
    print("Copying output...")
    subprocess.run(comando_copia, shell=True)
    print("Done")

#Copia la configurazione di output  nella configurazione di input
def copia_config(filename):
    comando_copia = "cp ../OUTPUT/CONFIG/config.spin "+filename
    print("Copying configuration...")
    subprocess.run(comando_copia, shell=True)
    print("Done")

#Modifica il tipo della simulazione (campo esterno nullo)
def modifica_sim_type(filename, value):
    with open(filename, 'r') as file:
        linee = file.readlines()

    with open(filename, 'w') as file:
        for linea in linee:
            if linea.startswith("SIMULATION_TYPE"):
                # Modifica il valore numerico associato all'intestazione "SIMULATION_TYPE"
                nuova_linea = f"SIMULATION_TYPE        {value} 1.0 0.0\n"
                file.write(nuova_linea)
            else:
                # Scrivi le altre righe del file senza modificarle
                file.write(linea)

if __name__ == "__main__":

    inputfile = "../INPUT/input.dat"
    configfile = "../INPUT/CONFIG/config.spin"
    methods = ["Metro", "Gibbs"]
    temps=np.linspace(0.5, 2.0, 40)
    temps_rev=temps[::-1]
    #print(temps_rev)
    sim_type=2
    for method in methods:
        modifica_restart(inputfile, 0)
        modifica_sim_type(inputfile, sim_type)
        for T in temps_rev:
            print("Simulating with "+method+" T=", format(T, '.3f'), "...")
            modifica_valore_input(inputfile, T)
            avvia_simulatore()
            modifica_restart(inputfile, value=1)
            copia_config(configfile)
            copia_output(T, method)
            print("Done")
        sim_type+=1
    
