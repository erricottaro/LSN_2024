import subprocess
import numpy as np

#avvia la simulazione NVE per esercizio 4
def avvia_simulatore():
    # Comando Bash per eseguire il programma simulator.exe
    comando_bash = "./simulator.exe"
    
    # Esegue il comando Bash
    print("Running simulation...")
    subprocess.run(comando_bash, shell=True)
    print("Done")

#Modifica la densità
def modifica_densita_input(filename, nuovo_valore):
    with open(filename, 'r') as file:
        linee = file.readlines()

    with open(filename, 'w') as file:
        for linea in linee:
            if linea.startswith("RHO"):
                # Modifica il valore numerico associato all'intestazione "T"
                nuova_linea = f"RHO                    {nuovo_valore}\n"
                file.write(nuova_linea)
            else:
                # Scrivi le altre righe del file senza modificarle
                file.write(linea)

#Modifica r_cut
def modifica_r_cut_input(filename, nuovo_valore):
    with open(filename, 'r') as file:
        linee = file.readlines()

    with open(filename, 'w') as file:
        for linea in linee:
            if linea.startswith("R_CUT"):
                # Modifica il valore numerico associato all'intestazione "T"
                nuova_linea = f"R_CUT                  {nuovo_valore}\n"
                file.write(nuova_linea)
            else:
                # Scrivi le altre righe del file senza modificarle
                file.write(linea)

#Modifica delta nel file di input
def modifica_delta_input(filename, nuovo_valore):
    with open(filename, 'r') as file:
        linee = file.readlines()

    with open(filename, 'w') as file:
        for linea in linee:
            if linea.startswith("DELTA"):
                # Modifica il valore numerico associato all'intestazione "T"
                nuova_linea = f"DELTA                  {nuovo_valore}\n"
                file.write(nuova_linea)
            else:
                # Scrivi le altre righe del file senza modificarle
                file.write(linea)

#Modifica la temperatura nel file di input
def modifica_temp_input(filename, nuovo_valore):
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

#clean output
def clean_output():
    comando_rimuovi = "rm ../OUTPUT/*.dat"
    subprocess.run(comando_rimuovi, shell=True)

#Copia i file di output in apposite cartelle ({fase})
def copia_output(phase):
    dir_name = f"../../{phase}"
    comando_dir = "mkdir -p "+dir_name
    comando_copia = "cp ../OUTPUT/*.dat "+dir_name
    print("Building directory...")
    subprocess.run(comando_dir, shell=True)
    print("Done")
    print("Copying output...")
    subprocess.run(comando_copia, shell=True)
    print("Done")

#Copia il file di input nella cartella di output
def copia_input(phase):
    dir_name = f"../../{phase}"
    comando_dir = "mkdir -p "+dir_name
    comando_copia = "cp ../OUTPUT/*.dat "+dir_name
    print("Building directory...")
    subprocess.run(comando_dir, shell=True)
    print("Done")
    comando_copia = f"cp ../INPUT/input.dat {dir_name}/input.dat"
    print("Copying input...")
    subprocess.run(comando_copia, shell=True)
    print("Done")

if __name__ == "__main__":

    inputfile = "../INPUT/input.dat"
    phases = ["Solid", "Liquid", "Gas"]
    rhos = [1.1, 0.8, 0.05]
    #thermalization temperatures
    temps = [1.55, 1.98, 0.9]
    r_cuts = [2.2, 2.5, 5.0]
    deltas = [0.001, 0.0005, 0.0005]
    clean_output()
    for i in range(3):
        rho = rhos[i]
        temp = temps[i]
        r_cut = r_cuts[i]
        phase = phases[i]
        delta = deltas[i]
        print(phase)
        modifica_densita_input(inputfile, rho)
        modifica_temp_input(inputfile, temp)
        modifica_r_cut_input(inputfile, r_cut)
        modifica_delta_input(inputfile, delta)
        copia_input(phase)
        #avvia simulazione
        avvia_simulatore()
        copia_output(phase)