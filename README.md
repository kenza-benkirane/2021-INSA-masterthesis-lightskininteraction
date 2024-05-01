# Light-Skin Interaction Simulation for PhotoPlethysmoGraphiy (PPG)

This project is my masters' thesis within the INSA Lyon Engineering Degree, held at Infineon Technologies AG, as part of the UNITECH International Engineer scheme. It was held from March to August 2021 in Munich, Germany.

More info on the company and there effort towards smart-wearables is available [here]([url](https://www.infineon.com/cms/en/applications/healthcare-and-lifestyle/wearables/)).

# Introduction and context
Wearable’s devices are questioning how people track their health and fitness. Recent studies have shown that, while elderly and sick people used to dish stop using eHealth wearables after 30 to 60 days, the global pandemic deeply changed their lifestyles and the way they track their health, and this is most likely to last. 
These main changes are particularly real when it comes to cardiac problems: prevention plays a key role in cardiac monitoring. A high HR in a healthy population is linked to a higher risk of coronary artery disease and can be a predictor of Heart Failure (HF). When wrongly recovered, HF can increase adverse cardiovascular events. The key problem about HR and HF, is that it arises in healthy patients, hence its’ nickname “silent killer”, and is also one of the first causes of death, especially in Europe. In the US, it accounts for 1 in 4 deaths. 
To monitor health quality, a consistent follow-up is deeply needed especially for conditions such as myocardial infarction or more complex arrhythmias. Now, persons with risks are followed through ECG and their cardiologist. It’s particularly interesting for patients with risk factors such as high BP, chest pain, shortness of breath,etc. The main problem remains that ECG is useful during a heart attack, which cannot be predicted and needs an appointment with the physician. Plus, the chances of surviving a heart attack are better the sooner the emergency treatment arises. Moreover, two-thirds of all heart attacks go undetected by ECG. 
That’s when PPG comes: PhotoPlethysmoGraphic (PPG) method is a great solution for personalized e-health: It allows following the oxygen rate, blood pulse and multiple other cardiac parameters by just putting a wristband or putting the device around one’s finger. It stands at the future ECG replacement. Moreover, it’s super user-friendly: it’s simple, has a great cost-benefit ratio, is easily accessible, non-invasive, and the signal acquisition is easy to read. This tool is especially interesting when it comes to detecting cardiac arrest or haemodynamic shock: it could save lives. The greatest advantage is that PPG monitoring is just one click away, which would alert and support the identification of the heart failure problem. Also, some smart wristbands offer an emergency call option: These features make the PPG more convenient and adapted to reactive support to the patient. 
This project takes part in a long range of research projects within the Innovation Department at Infineon Technology. It takes great inspiration in previous results and was inspired by the Non-Destructive Testing module taken at my engineering school, INSA Lyon. 

# How it works
The projects aims at simulating the light-skin interaction between a PPG and human skin: the simulation focuses on the skin from a layered material point of view and analyses the interaction and absorption-reflection dynamics with lights, depending on the wavelength and intensity of the LEDs.
This code was created for simulating Diffuse Reflectance Spectroscopy(DRS) of human skin tissue via Monte Carlo for Multi-Layered media(MCML).You may need to refer to this paper to understand the principle: https://omlc.org/software/mc/mcml/index.html and https://omlc.org/software/mc/mcml/MCman.pdf

<img width="760" alt="Screenshot 2023-10-21 at 11 01 42 PM" src="https://github.com/kenzabenki/INSA-masterthesis-lightskininteraction/assets/52839072/36c788a1-86df-490e-ab14-8089152884e7">


# Usage
These file must be in same path with **run.py**</br>
+ wavelength.csv: all wavelength data for this program
+ model_input.txt: some parameters of the simulated tissue model
+ mua_water.csv: the optical absorption properties (μ<sub>a</sub> cm<sup>-1</sup>) of water
+ mua_oxy.csv: the optical absorption properties (μ<sub>a</sub> cm<sup>-1</sup>) of oxygenated whole blood
+ mua_deoxy.csv: the optical absorption properties (μ<sub>a</sub> cm<sup>-1</sup>) of deoxygenated whole blood
+ mua_melanin.csv: the optical absorption properties (μ<sub>a</sub> cm<sup>-1</sup>) ofinterior of typical cutaneous melanosome</br>

And these optical absorption properties referred to https://omlc.org/software/mc/mcxyz/index.html
- - -
This code is a preset 8-multicore processing. If you have enough CPU cores, you must increase/revise the code.</br>
You must revise [line731](https://github.com/GarrettTW/MCML_simulate-spectroscopy/blob/33d8c457c14ce4164e525d4fda282cfcbaf2abc0/run.py#L731) the **8** to **number of your CPU cores**

```python
if cpu_number>=1 and cpu_number<=8 and cpu_number<=cpu_count(): 
```

You must increase the following code below [line791](https://github.com/GarrettTW/MCML_simulate-spectroscopy/blob/33d8c457c14ce4164e525d4fda282cfcbaf2abc0/run.py#L791) (if you want to use 9 cores for mutiprocessing)

```python
if c.get(8):
    boundary = c[8]
    q9 = Queue()
    m9 = mp.Process(target=job, args=(q9,model,N,boundary))
    m9.start()                        
```
and increase the following code below [line824](https://github.com/GarrettTW/MCML_simulate-spectroscopy/blob/33d8c457c14ce4164e525d4fda282cfcbaf2abc0/run.py#L824)

```python
if c.get(8):
    m9.join()
    R.append(q9.get())                     
```
- - -
Run the `run.py`

Enter your requirements as prompted

`save the Reflection Spectrum?(y/n):`You can type`y`or`n`to choose whether to store the simulated spectrum on your hard drive.
    
`the filename?(do NOT add extension):`If you type`y`, you need to type the file name (you do **not** need to enter the **filename extension**)
    
`How many photons for simulation?(1000 photons spend about 4 mins)`more photons will result in better spectrum (but you will spend more time)

`How many multicore operations to use?`Please consider the number of your CPU cores

# Ressources
This code was created for simulating Diffuse Reflectance Spectroscopy(DRS) of human skin tissue via Monte Carlo for Multi-Layered media(MCML)
You may need to refer to this paper to understand the principle.

In all simulations, 10<sup>4</sup> photon packets at each wavelength were used to enter the model from above (z<0) at the air at stratum corneum boundary. At this point a set of probabilities determine if the photon packets are reflected or enter the tissue. When a packet has entered the tissue, it can be partially or fully absorbed by event bins uniformly distributed throughout the tissue.With a starting photon weight W of 1, every time a photon packet interacts with a bin it loses part of its weight and then gets scattered in a direction determined by the anisotropy factor and scattering coefficient. At the end of the simulation, all diffusely scattered photons locating at the incident side (z<0) were added up to give a diffuse reflectance spectrum. </br>
One of the simulation results is presented in **below figure** with filled circles enclosed in error bars estimated from 10 simulations.
For comparison, a diffuse reflectance spectrum measured on the inside ofthe right forearm of an adult with our DRS apparatus is also presented (red solid curve). A good agreement with the simulation profile was obtained, suggesting that the model used captured the essential elements of human skin tissue.
![Imgur](https://i.imgur.com/cHXjQje.jpg "Monte-Carlo multilayer simulation (solid circles with error bars) and measured (red solid curve) diffuse reflectance spectrum of human skin tissue.")
