import subprocess
import numpy as np

#avvia la simulazione a temperatura T
def avvia_simulatore():
    # Comando Bash per eseguire il programma es08.2.x
    comando_bash = "./es08.2.x"
    
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

#Copia i file di output in apposite cartelle
def copia_output(temp):
    dir_name = f"T={format(temp, '.3f')}"
    comando_dir = "mkdir -p "+dir_name
    comando_copia = "cp OUTPUT/annealing_energy.dat OUTPUT/annealing_accep.dat OUTPUT/annealing_output.dat "+dir_name
    print("Building directory...")
    subprocess.run(comando_dir, shell=True)
    print("Done")
    print("Copying output...")
    subprocess.run(comando_copia, shell=True)
    print("Done")

#copia mu e sigma nel system_input
def copia_params(filename, mu, sigma):
    with open(filename, 'r') as file:
        linee = file.readlines()

    with open(filename, 'w') as file:
        for linea in linee:
            if linea.startswith("MU"):
                # Modifica il valore numerico associato all'intestazione "MU"
                nuova_linea = f"MU                   {mu}\n"
                file.write(nuova_linea)
            elif linea.startswith("SIGMA"):
                # Modifica il valore numerico associato all'intestazione "SIGMA"
                nuova_linea = f"SIGMA               {sigma}\n"
                file.write(nuova_linea)
            else:
                # Scrivi le altre righe del file senza modificarle
                file.write(linea)


if __name__ == "__main__":
    ann_input = "INPUT/annealing_input.dat"
    system_input = "INPUT/system_input.dat"
    #reset parametri a 0.8 e 0.5
    copia_params(system_input, 0.8, 0.5)
    for i in range(21):
        temp = (i+2.0)*np.exp(-0.5*i)
        print("=======================")
        print("T=", format(temp, '.3f'))
        modifica_valore_input(ann_input, temp)
        avvia_simulatore()
        copia_output(temp)
        new_mu = np.loadtxt("final_parameters.dat", usecols=1, skiprows=0, max_rows=1, unpack=True)
        new_sigma = np.loadtxt("final_parameters.dat", usecols=1, skiprows=1, max_rows=1, unpack=True)
        print(new_mu, new_sigma)
        copia_params(system_input, new_mu, new_sigma)