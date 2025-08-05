# microtonal-controller
Firmware for microtonal keyboard controller

This is the firmware for the Desiderata Systems Seven Limit Mosaichord just-intonation keyboard controller and any other projects based on the same architecture.

https://desideratasystems.com

![Seven Limit](doc/seven_limit.jpeg?raw=true "Seven Limit controller")

The design currently uses a Teensy 4.0 or 4.1 microcontroller.  The recommended way to edit the source code and load new firmware is through the Arduino IDE with Teensyduino extensions.  (I'm using 2.3.2 of the Arduino IDE.)

My plan is for the firmware source code to remain open source.  The PCB layout will probably remain proprietary, but this repo has schematics.

Building alternative control surfaces should be relatively straightforward with simple parts.  The control surface bus is meant to be extensible, so you can daisy chain multiple control surfaces off of one controller module if you want (at some cost in the form of electrical noise and a slower scan rate).

The user manual covers how to compile the source code and load it onto the microcontroller.

As of version 1.1, we've added a local synthesizer that uses the Mosaichord's built-in DAC, so it's now possible to use it as a stand-alone instrument rather than just a controller.  (It uses a fairly typical subtractive architecture with basic waveforms into a per-voice resonant low-pass filter.)

The Mosaichord project improves upon an earlier prototype described in great detail here: http://jsnow.bootlegether.net/jik/keyboard.html
