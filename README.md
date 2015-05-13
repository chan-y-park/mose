# mose

MOduli Space Explorer: Python program to study moduli spaces of elliptic fibrations.

## To Do List

## How to run

### From Python interpreter
* Generating a K-wall network data
  1. In the root directory that contains ```mose``` directory, start the python interpreter.
  1. ```>>> import mose```
  1. Load a configuration file.
    1. Load it interactively.
      * ```>>> config = mose.load_config()```
      * A file dialog window opens, select a configuration file to use.
    1. Or give the file name as an argument.
      * ```>>> config = mose.load_config('default.ini')```.
    1. The returned value is ```MoseConfig``` that contains all the configuration information.
  1. Run ```analysis```, which returns ```data```.
    1. To get a K-wall network at a fixed phase, run
      * ```>>> data = mose.analysis(config, phase=1.0)```
    1. To get multiple K-wall networks according to the configuration, run
      * ```>>> data = mose.analysis(config)```
* Plotting a K-wall network
  1. Run ```make_plots```.
    * ```>>> mose.make_plots(config, data)```
  1. Click an object in each figure to display its label.
  1. Press ```d``` to delete all the displayed labels.
* Saving the data
  1. ```>>> mose.save(config, data)```
  1. Then you get a directory dialog window to select a directory to save ```config.ini``` and ```data.mose```.
* Loading the data
  1. Load it interactively.
    * ```>>> config, data = mose.load()```
    * Then you first get a directory dialog window to select a directory that contains the configuration & the data file. Select a directory and it returns ```(config, data)```.
  1. Or give the directory name as an argument.
    * ```>>> config, data = mose.load('data/name')```.

### From the Linux shell
* In the root directory, there is an executable Python script ```main.py```. Running it without any option prints the help screen.
* For example, if you want to run ```mose``` with
  * ```default.ini``` configuration file, 
  * for a fixed phase at 1.0, 
  * producing its plot,
  * and saving the resulting data at ```data/yyyy-mm-dd-tt-mm```,
```
> ./main.py --config-file=default.ini --phase=1.0 --show-plot --save-data
```
* If you want to load a data and show the plot, 
```
> ./main.py --load-data --show-plot
```
  1. Then you first get a directory dialog window to select a directory that contains the configuration & the data file. Select a directory and ```mose``` plots the K-wall network.
