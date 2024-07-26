# flightpath_reconstruction_using_ukf
In this project, flightpath reconstruction was carried out using unscented kalman filter. Cessna 172 was manually controlled in simulated environment then flight data was extracted from logged Flightgear data to perform FPC algorithm.

"protocol.xml" was used according to "https://wiki.flightgear.org/Generic_protocol". Running the Flightgear from command prompt with the "fgfs --config=protocol.xml" will log files into "example.csv" located (for Windows) in "\AppData\Roaming\flightgear.org\Export".

Here is resulting reconstructed flightpath:
![3D_plot](https://github.com/user-attachments/assets/0d7dfed5-f4e6-4c01-9274-2148a15eb181)

Unscented kalman filter was based on methods provided by following references:

Wan, Eric, and Rudolph Van Der Merwe. The Unscented Kalman Filter. 2004.

Chapter 7: Recursive Parameter Estimation 
"Flight Vehicle System Identification - A Time Domain Methodology"
Second Edition
by Ravindra V. Jategaonkar
published by AIAA, Reston, VA 20191, USAA

Otávio, Bruno & Teixeira, Soares & Torres, Leonardo & Henriques, 
Paulo & Oliveira, Iscold & Aguirre, Luis. (2005). 
Flight path reconstruction using the unscented Kalman filter algorithm. 
