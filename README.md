# microtonal-controller
Firmware for microtonal keyboard controller

This is the firmware for the JI7L4 just-intonation keyboard controller and any other projects based on the same architecture.

The design currently uses a Teensy 4.0 or 4.1 microcontroller.  The recommended way to edit the source code and load new firmware is through the Arduino IDE with Teensyduino extensions.  (I'm using 2.1.0 of the Arduino IDE.)

The project is in an early state and not yet available generally yet, but my plan is for the firmware source code as well as the schematic to be open source.  The PCB designs will probably remain proprietary.

Building alternative control surfaces should be relatively straightforward with simple parts.  The control surface bus is meant to be extensible, so you can daisy chain multiple control surfaces off of one controller module if you want (at some cost in the form of electrical noise and a slower scan rate).
