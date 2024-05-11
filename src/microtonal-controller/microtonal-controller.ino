/**@file microtonal-controller.ino */
/*
Copyright 2023 Jim Snow

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#define hwversion 3

#include <MIDI.h>

#include <ILI9341_t3.h>
#include <font_Arial.h>
#include <font_ArialBold.h>

#include <ADC.h>
#include <WS2812Serial.h>
#include <FlexCAN_T4.h>


/* Pins */

#if (hwversion <= 1)
#define ledPin 1
#define adc1Pin 15
#define adc2Pin 16
#define adc3Pin 17
#define adc4Pin 18
#define shiftRegisterClockPin 5
#define shiftRegisterOutPin 4
#define midiOutPin 20
#define midiInPin 0
#define sckPin 13
#define sdiPin 11
#define sdoPin 12
#define screenCSPin 10
#define backlightPin 3
#define screenDCPin 9
#define touchCSPin 8
#endif

#if (hwversion == 2)
#define mutePin 0
#define touchIrqPin 1
#define i2sOutPin 2
#define i2sLRClkPin 3
#define i2sBClkPin 4
#define shiftRegisterClockPin 5
#define touchCSPin 6
#define shiftRegisterOutPin 7
#define ledPin 8
#define screenDCPin 9
#define screenCSPin 10
#define sdiPin 11
#define sdoPin 12
#define sckPin 13
#define adc4Pin 15
#define adc3Pin 16
#define adc2Pin 17
#define adc1Pin 18
#define backlightPin 19
#define midiOutPin 20
#define midiInPin 21
#endif

#if (hwversion < 3)
#define pullupPin1 255
#define pullupPin2 255
#define pullupPin3 255
#define pullupPin4 255
#endif

#if (hwversion == 3)
#define pullupPin1 0
#define pullupPin2 1
#define i2sOutPin 2
#define i2sLRClkPin 3
#define i2sBClkPin 4
#define shiftRegisterClockPin 5
#define pullupPin3 6
#define shiftRegisterOutPin 7
#define ledPin 8
#define screenDCPin 9
#define pullupPin4 10
#define sdiPin 11
#define sdoPin 12
#define sckPin 13
#define adc4Pin 15
#define adc3Pin 16
#define adc2Pin 17
#define adc1Pin 18
#define backlightPin 19
#define midiOutPin 20
#define midiInPin 21
#define canTxPin 22
#define canRxPin 23

#define screenCSPin 36
#endif


uint32_t debugFlags = 0;

enum dFlags {
  adcCalibrationDebug = 0,
};

#define dbg(flag) ((debugFlags & (1 << flag)) != 0)
#define dbgSet(flag) debugFlags = debugFlags | (1 << flag)
#define dbgClear(flag) debugFlags = debugFlags & (~(1 << flag))

/* Serial */

void serialSetup() {
  Serial.begin(115200);
}


/* LEDs */

const int numLeds = 6;
byte drawingMemory[numLeds * 3];
DMAMEM byte displayMemory[numLeds * 12];
WS2812Serial leds(numLeds, displayMemory, drawingMemory, ledPin, WS2812_GRB);

#define RED    0xFF0000
#define GREEN  0x00FF00
#define BLUE   0x0000FF
#define YELLOW 0xFFFF00
#define PINK   0xFF1088
#define ORANGE 0xE05800
#define WHITE  0xFFFFFF
#define BLACK  0x000000
#define DIM    0x000100

/*
#define RED    0x160000
#define GREEN  0x001600
#define BLUE   0x000016
#define YELLOW 0x101400
#define PINK   0x120009
#define ORANGE 0x100400
#define WHITE  0x101010
*/

void ledSetup() {
  leds.begin();

  int microsec = 15000 / leds.numPixels();
  
  colorWipe(RED, microsec);
  colorWipe(GREEN, microsec);
  colorWipe(BLUE, microsec);
  colorWipe(YELLOW, microsec);
  colorWipe(PINK, microsec);
  colorWipe(ORANGE, microsec);
  colorWipe(WHITE, microsec);
  colorWipe(BLACK, microsec);
  leds.setPixel(0, DIM);
  leds.show();
}

int colorDiv (int in, int div) {
  int r = (in & 0xff0000) >> 16;
  int g = (in & 0x00ff00) >> 8;
  int b = in & 0x0000ff;

  return (r/div << 16) | (g/div << 8) | (b/div);
}

void colorWipe(int color, int wait) {
  for (int i=0; i < leds.numPixels(); i++) {
    leds.setPixel(i, colorDiv(color, 64));
    leds.show();
    delayMicroseconds(wait);
  }
}

/* AUDIO */

#include <Audio.h>
#include <Wire.h>
#include <SPI.h>
#include <SD.h>
#include <SerialFlash.h>

// GUItool: begin automatically generated code
AudioSynthWaveformSine   sine1;          //xy=323.3333396911621,295.3333330154419
AudioOutputI2S2          i2s2_1;         //xy=597.3333282470703,333.3333282470703
AudioConnection          patchCord1(sine1, 0, i2s2_1, 0);
AudioConnection          patchCord2(sine1, 0, i2s2_1, 1);
AudioControlSGTL5000     sgtl5000_1;     //xy=445.3333282470703,400.3333282470703
// GUItool: end automatically generated code

/* A440 test tone */

void audioSetup() {
  AudioMemory(20);
  AudioNoInterrupts();
  sine1.frequency(440.0);
  sine1.amplitude(0.0);
  AudioInterrupts();
  //noise1.amplitude(1.0);
  //pinMode(mutePin, OUTPUT);
  //digitalWrite(mutePin, HIGH); /* unmute */
  Serial.println("audio setup complete");
}

void setPitchReference(double freq) {
  AudioNoInterrupts();
  sine1.frequency(freq);
  sine1.amplitude(0.1);
  AudioInterrupts();
}

/* ADCs */

#define adcChannels 4
const int adcPins[adcChannels] = {adc1Pin, adc2Pin, adc3Pin, adc4Pin};
const int adcPullupPins[adcChannels] {pullupPin1, pullupPin2, pullupPin3, pullupPin4};

ADC *adc = new ADC();

void adcSetup() {

  pinMode(LED_BUILTIN, OUTPUT);
  for (int i=0; i<adcChannels; i++) {
    int adcPin = adcPins[i];
    int pullupPin = adcPullupPins[i];
    pinMode(adcPin, INPUT_DISABLE);
    if (pullupPin != 255) {
      pinMode(pullupPin, OUTPUT);
      digitalWrite(pullupPin, HIGH);
    }
  }

  auto convSpeed = ADC_CONVERSION_SPEED::HIGH_SPEED;
  auto sampleSpeed = ADC_SAMPLING_SPEED::HIGH_SPEED;
  int resolution = 12;
  int averaging = 4;

  adc->adc0->setAveraging(averaging);
  adc->adc0->setResolution(resolution);
  adc->adc0->setConversionSpeed(convSpeed);
  adc->adc0->setSamplingSpeed(sampleSpeed);

  adc->adc1->setAveraging(averaging);
  adc->adc1->setResolution(resolution);
  adc->adc1->setConversionSpeed(convSpeed);
  adc->adc1->setSamplingSpeed(sampleSpeed);
}

int otherChannel(int thisChannel) {
  switch (thisChannel) {
    case (0): return 2;
    case (1): return 3;
    case (2): return 0;
    case (3): return 1;
  }
  return 0;
}

/* 
 * Determine how many microseconds to pause before reading ADCs,
 * to give cicuit time to settle after advancing the shift register.
 * We look at the immediate previous 4 values read, and the values
 * for the current 4 as of the last update.
 * If the difference is large on any channel, then we wait longer.
 */
int getADCDelay(int *prev, int *curr) {
  int maxDelay = 12;
  for (int channel = 0; channel < 4; channel++) {
    int delta = prev[channel] - curr[channel];
    if (delta < 0) {
      delta = -delta;
    }

    int delay = delta / 100;
    if (delay > maxDelay) {
      maxDelay = delay;
    }
  }
  return maxDelay;
}

float adcScale = 1.0f/4095.0f;

/*
 * determine resistance that causes a voltage sag on ADC inputs pulled up
 * by 3.3 volts with a 3.3k or 1k resistor, 200 ohm series resistor
 */
float valueToResistance(int value) {
  float v = value * adcScale;
#if (hwversion < 3)
  float pullup = 3300.0f;
#else
  float pullup = 1000.0f + 20.0f; /* default pin output impedance is probably around 20 ohms */
#endif
  float series = 206.0f; /* 200 ohm resister, about 6 ohms more for shift register */

  /* avoid divide by zero */
  if (1.0f - v < 0.01f) {
    v = 1.0f - 0.01f;
  }

  float r = ((v * pullup) / (1.0f - v) - series);

  if (r < 1.0) {
    r = 1.0;
  }
  return r;
}

/*
 * Measure electrical resistance between channels across the FSR.  Low
 * resistance can throw off ADC readings, so we have to compensate.
 * The resistances change depending on what keys are being pressed, so
 * we need to re-calibrate on every keyboard scan.
 *
 * This should be called before the shift registers have selected bit 0,
 * when the addc channel readings are only affected by driving the pullup
 * pins high or low.
 */
void calibrateADCs(bool verbose, float cal[adcChannels][adcChannels]) {
  for (int i = 0; i < adcChannels; i++) {
    for (int j = 0; j < i; j++) {
      if (i==j) {
        cal[i][j] = 1000000.0f;
        continue;
      }

      for (int channel = 0; channel < adcChannels; channel++) {
        int pin = adcPullupPins[channel];
        if (channel == i) {
          pinMode(pin, OUTPUT);
          digitalWrite(pin, HIGH);
        } else if (channel == j) {
          pinMode(pin, OUTPUT);
          digitalWrite(pin, LOW);
        } else {
          pinMode(pin, INPUT_DISABLE);
        }
      }

      delayMicroseconds(50);

      int v1,v2;

      readADCs(i, v1, j, v2);

      float iadc = ((float)v1) * adcScale;
      float jadc = ((float)v2) * adcScale;

      float iv = 1.0f - (((1.0f - iadc) * 1220.0f) / 1020.0f);
      float jv = (jadc * 1220.0f) / 1020.0f;
      float rijTotal = (2440.0f / (1.0f - (iv -jv))) - 2440.0f;

      //int k, l;
      //getKl(i, j, k, l);

      //readADCs(k, v1, l, v2);

      //float kadc = ((float)v1) * adcScale;
      //float ladc = ((float)v2) * adcScale;

      if (verbose) Serial.println("calibration i:" + String(i) + " j:" + String(j) + " iv:" + String((1.0f - iv) * 100.0f) + " jv:" + String(jv * 100.0f) + " rijTotal:" + String(rijTotal));

      //if (verbose) {
      //  Serial.println("calibration i:" + String(i) + " j:" + String(j) + " iv:" + String((1.0f - iv) * 100.0f) + " jv:" + String(jv * 100.0f) + " rijTotal:" + String(rijTotal) +
      //    " k:" + String(k) + "(" + String(kadc * 100.0f) + ") l:" + String(l) + "(" + String(ladc * 100.0f) + ")");
      //}
      /* store the reciprocal to we don't have to do division later */
      cal[i][j] = cal[j][i] = 1.0f / rijTotal;
    }
  }

  for (int i = 0; i < adcChannels; i++) {
    int pin = adcPullupPins[i];
    pinMode(pin, OUTPUT);
    digitalWrite(pin, HIGH); /* must come after setting mode to output, otherwise ignored */
  }

  delayMicroseconds(40);
}

/*
 * Convert a triangle of resistors into an equivalent set of three resistors
 * that meet in a junction.  Uses reciprocal resistance.
 */
void deltaYConv(float ab, float bc, float ac, float& ra, float& rb, float &rc) {
  ra = ab + ac + (ab * ac) / bc;
  rb = ab + bc + (ab * bc) / ac;
  rc = ac + bc + (bc * ac) / ab;
}

void rsScale(int i, int j, float scale, float weight, float rs[adcChannels][adcChannels]) {
  float orig = rs[i][j];
  float val = (orig * (1.0f - weight)) + (orig * scale * weight);

  rs[i][j] = val;
  rs[j][i] = val;
}

/*
 * Compute resistance across two nodes i and j in a fully-connected network of four nodes.
 * Uses reciprocal resistances.
 *
 * This isn't solvable by simplifying parallel and series resistors, so we use delta-Y
 * conversion to transform triangle ilk into three resistors meeting at a new point m.
 *
 * Then we can solve it as series and parallel resistances.
 */
float rij(int i, int j, int k, int l, float rs[adcChannels][adcChannels]) {
  float ij = rs[i][j];
  float ik = rs[i][k];
  float il = rs[i][l];
  float jk = rs[j][k];
  float jl = rs[j][l];
  float kl = rs[k][l];

  float im, lm, km;
  deltaYConv(il, kl,ik, im, lm, km);

  float jm1 = 1.0f / ((1.0f / jl) + (1.0f / lm));
  float jm2 = 1.0f / ((1.0f / jk) + (1.0f / km));

  float jm = jm1 + jm2;

  float ij2 = 1.0f / ((1.0f / im) + (1.0f / jm));

  return ij + ij2;
}

/*
 * The calibration routine finds the resistance between each channel to the other.
 * That's not exactly what we want, though -- we need to know what the resistance between
 * any two channels would be if we could ignore the electrical path through the other two
 * channels.
 *
 * We figure this out by making an initial guess for the values of all the phantom resistors,
 * calculating what the resistance should be from one channel to the next through all available
 * paths, and then adjusting the values of the resistors according to the discrepency between
 * the two. Hopefully we converge on a solution.
 *
 * There may be more than one solution, in which case our guess will probably be somewhat
 * off.
 *
 * Resistance values are stored as 1/r to avoid division.
 */
void refineCalibration (const float in[adcChannels][adcChannels], float out[adcChannels][adcChannels], int iterations) {
  float avg;
  float sum = 0.0f;
  for (int i = 1; i < adcChannels; i++) {
    for (int j = 0; j < i; j++) {
       sum += in[i][j];
    }
  }

  avg = sum / 6.0f;

  /*
   * Initial guess: adjacent channels have higher throughput than non-adjacent channels, and input matrix
   * resistance is a lot less than it should be.
   *  */
  for (int i = 0; i < adcChannels; i++) {
    for (int j = 0; j < adcChannels; j++) {
      float val;
      float scale = 0.5f;
      if (i==j) {
        val = 1000000.0f;
      } else if ((i + j) % 2 == 0) {
        val = avg * 0.5f * scale;
      } else {
        val = avg * 1.25f * scale;
      }
      out[i][j] = val;
    }
  }

  while (iterations-- > 0) {
    for (int i = 0; i < adcChannels; i++) {
      for (int j = 0; j < adcChannels; j++) {
        if (i==j) {
          continue;
        }

        int k, l;
        getKl(i, j, k, l);

        float ijGuess = rij(i, j, k, l, out);
        float ijActual = in[i][j];

        float scale = ijActual / ijGuess;

        rsScale(i, j, scale, 1.0f, out);
        rsScale(i, k, scale, 0.5f, out);
        rsScale(i, l, scale, 0.5f, out);
        rsScale(j, k, scale, 0.5f, out);
        rsScale(j, l, scale, 0.5f, out);
        /* leave k,l as it is */

        if (dbg(adcCalibrationDebug)) {
          Serial.println("refineCalibration " + String(i) + " " + String(j) + " " + String(1.0f/ijGuess) + " " + String(1.0f/ijActual));
        }
      }
    }
  }
}

/* given a number from 0-3, populate other arguments with remaining digits in arbitrary order */
void getJkl(const int i, int &j, int &k, int &l) {
  for (j = 0; j < adcChannels; j++) {
    if (j != i) {
      break;
    }
  }
  for (k = 0; k < adcChannels; k++) {
    if (k != i && k != j) {
      break;
    }
  }
  for (l = 0; l < adcChannels; l++) {
    if (l != i && l != j && l !=k) {
      break;
    }
  }

  if (i + j + k + l != (adcChannels * (adcChannels-1)) / 2) {
    Serial.println("getJkl logic error");
  }
}

void getKl(const int i, const int j, int &k, int &l) {
  for (k = 0; k < adcChannels; k++) {
    if (k != i && k != j) {
      break;
    }
  }

  for (l = 0; l < adcChannels; l++) {
    if (l != i && l != j && l !=k) {
      break;
    }
  }

  if (i + j + k + l != (adcChannels * (adcChannels-1)) / 2) {
    Serial.println("getKl logic error");
  }
}

/*
 * Compute sum of voltages 5, each through its own respective resistor.
 * To avoid divisions, resistor values are input as reciprocals (1.0/r).
 */
inline float avg5(float v1, float r1, float v2, float r2, float v3, float r3, float v4, float r4, float v5, float r5) {
  return (v1*r1 + v2*r2 + v3*r3 + v4*r4 + v5*r5) / (r1 + r2 + r3 + r4 + r5);
}

inline float avg2 (float v1, float r1, float v2, float r2) {
  return (v1*r1 + v2*r2) / (r1 + r2);
}

/* 
 * Like avg5, but solve for r1 if we know average voltage already (v).
 * Resistances again input as reciprocals, and output is likewise.
 * derivation from avg5:
 *
 * (v1*r1 + v2*r2 + v3*r3 + v4*r4 + v5*r5) / (r1 + r2 + r3 + r4 + r5) = v
 * (v1*r1 + v2*r2 + v3*r3 + v4*r4 + v5*r5) = v*r1 + v*(r2 + r3 + r4 + r5)
 * v1*r1 + v2*r2 + v3*r3 + v4*r4 + v5*r5 - (v * (r2 + r3 + r4 + r5))  = v*r1
 * v2*r2 + v3*r3 + v4*r4 + v5*r5 - (v * (r2 + r3 + r4 + r5)) = v*r1 - v1*r1
 * v2*r2 + v3*r3 + v4*r4 + v5*r5 - (v * (r2 + r3 + r4 + r5)) = r1 (v-v1)
 * (v2*r2 + v3*r3 + v4*r4 + v5*r5 - (v * (r2 + r3 + r4 + r5)))/(v-v1) = r1
 */
inline float avg5r1(float v, float v1, float v2, float r2, float v3, float r3, float v4, float r4, float v5, float r5) {
  float r1 = (v2*r2 + v3*r3 + v4*r4 + v5*r5 - (v * (r2 + r3 + r4 + r5))) / (v - v1);

  if(r1 < 0.0) {
    return (1.0f / 1000000.0f);
  }

  /*
  if (isnan(r1)) {
    return 1.0f / 100000.0f;
  }*/
  return r1;
}

/*
 * Given a calibration matrix cal, voltages sampled by the ADC and some initial values for the resistances of four channels,
 * compute more accurate resistances.
 * This is only approximate, with "damping" as a tuning parameter.
 */
void applyCalibration(bool verbose, const float cal[adcChannels][adcChannels], const float vAdc[adcChannels], float r[adcChannels], const int iterations) {
  float v[adcChannels];
  float vNext[adcChannels];
  float rNext[adcChannels];
  float vOrig[adcChannels];
  //float rOrigInv[adcChannels];

  for (int i = 0; i < adcChannels; i++) {
    /* approximate first guess, compensate for 200 ohm resistor */
    v[i] = 1.0 - (( max((1.0 - vAdc[i]), 0.01)  * 1220.0) / 1020.0);
    vOrig[i] = v[i];
    //rOrigInv[i] = 1.0 / r[i];
  }

  if (verbose) {
    Serial.println(String(vAdc[0]) + " " + String(vAdc[1]) + " " + String(vAdc[2]) + " " + String(vAdc[3]) +
                   "|" + String(v[0]) + " " + String(v[1]) + " " + String(v[2]) + " " + String(v[3]) +
                   "|" + String(r[0]) + " " + String(r[1]) + " " + String(r[2]) + " " + String(r[3]));
  }

  float damping = 0.8;

  for (int iteration = 0; iteration < iterations; iteration++) {


    for (int i = 0; i < adcChannels; i++) {
      #if 1
      int j, k, l;
      getJkl(i, j, k, l);

      /* figure out what value of r1 would cause the voltage we're seeing at v[i] */
      float ri = 1.0f / avg5r1(vOrig[i], 0.0f, 1.0f, 1.0f/1220.0f, v[j], cal[i][j] * damping, v[k], cal[i][k] * damping, v[l], cal[i][l] * damping);
      rNext[i] = ri;
      /* now calculate what voltage v[0] ought to be if the other three channels weren't interfering */
      vNext[i] = avg2(1.0f, (1.0 / 1220.0f), 0.0, (1.0/ri) );

      //vNext[i] = avg5(1.0f, 1220.0f, v[j], 1.0/cal[i][j], v[k], 1.0/cal[i][k], v[l], 1.0/cal[i][l], 0.0, 1.0/rOrigInv[i]);

      #else
      vNext[i] = vOrig[i];
      for (int j = 0; j < adcChannels; j++) {
        if (i == j) {
          continue;
        }
        int k, l;
        getKl(i, j, k, l);

        float vdelta = v[j] - v[i];
        float strength = (cal[i][j]) / ((1.0f / 1220.0f) + (rOrigInv[i]) + (cal[i][j]) + (cal[i][k]) + (cal[i][l]) );

        vNext[i] -= vdelta * strength * damping;
      }
      #endif
    }

    /*
    for (int i = 0; i < adcChannels; i++) {
      rNext[i] = (1220.0f / (1.0f - v[i])) - 1220.0f;
    } */

    for (int i = 0; i < adcChannels; i++) {
      v[i] = min(max(vNext[i], 0.01f), 0.99f);
      r[i] = max(rNext[i], 1.0f);
    }
    if (verbose) {
      Serial.println("applyCalibration " + String(v[0]) + " " + String(v[1]) + " " + String(v[2]) + " " + String(v[3]) + " " + String(r[0]) + " " + String(r[1]) + " " + String(r[2]) + " " + String(r[3]));  
    }
  }
}

int readADC(int channel) {
  return adc->adc0->analogRead(adcPins[channel]);
}

void readADCs(int channel1, int &value1, int channel2, int &value2) {
  adc->adc0->startSingleRead(adcPins[channel1]);
  adc->adc1->startSingleRead(adcPins[channel2]);

  while(!adc->adc0->isComplete()) {};
  value1 = adc->adc0->readSingle();

  while(!adc->adc1->isComplete()) {};
  value2 = adc->adc1->readSingle();
}

void readADCs(int values[adcChannels]) {
  /* do two reads at a time in parallel */
  adc->adc0->startSingleRead(adcPins[0]);
  adc->adc1->startSingleRead(adcPins[1]);

  while(!adc->adc0->isComplete()) {};
  values[0] = adc->adc0->readSingle();
  while(!adc->adc1->isComplete()) {};
  values[1] = adc->adc1->readSingle();

  adc->adc0->startSingleRead(adcPins[2]);
  adc->adc1->startSingleRead(adcPins[3]);

  while(!adc->adc0->isComplete()) {};
  values[2] = adc->adc0->readSingle();
  while(!adc->adc1->isComplete()) {};
  values[3] = adc->adc1->readSingle();
}

void computeResistances(bool verbose, const int values[adcChannels], const float calibrationMatrix[adcChannels][adcChannels], float preCalibrationResistances[adcChannels], float resistances[adcChannels], int iterations) {
  float vAdc[adcChannels];

  for (int i = 0; i < adcChannels; i++) {
    vAdc[i] = values[i] * adcScale;
    float r = valueToResistance(values[i]);
    preCalibrationResistances[i] = r;
    resistances[i] = r;
  }

  applyCalibration(verbose, calibrationMatrix, vAdc, resistances, iterations);
}

float zeroPressureResistance = 5000.0f;
float maxPressureResistance = 600.0f;

/*
 * Force can be outside range of 0 (minimum force) to 1 (maximum force)
 */
inline float resistanceToForce(float r, float area = 1.0f) {
  return 1.0f - lerpNoClamp(max(area, 2.0f) / maxPressureResistance, 1.0f / r, area / zeroPressureResistance);
}

/* Shift Registers */

void shiftRegisterSetup() {
  pinMode(shiftRegisterOutPin, OUTPUT);
  pinMode(shiftRegisterClockPin, OUTPUT);
  shiftRegisterReset(0);
}

const int maxShiftRegisterBits = 8+32;

/* 
 * shift the clock until the whole register is cleared,
 * load a bit at the beginning, but don't update output
 * yet -- no bits are "visible" and curBit is implicitly -1
 */
void shiftRegisterReset(int curBit) {
  for (int i = curBit; i < maxShiftRegisterBits; i++) {
    shiftRegisterClock();
  }
  /* load a bit into the shift register */
  digitalWrite(shiftRegisterOutPin, HIGH);
  delayMicroseconds(2);
  shiftRegisterClock();
  digitalWrite(shiftRegisterOutPin, LOW);
}

/*
 * reset shift register and cycle the clock, so the output register shows bit 0 set
 * curBit is implicitly zero
 */
void shiftRegisterResetLoadBit(int curBit) {
  shiftRegisterReset(curBit);
  shiftRegisterClock();
}

void shiftRegisterClock() {
  digitalWrite(shiftRegisterClockPin, HIGH);
  delayMicroseconds(2);
  digitalWrite(shiftRegisterClockPin, LOW);
  delayMicroseconds(2);
}


/* Screen */

#define TFT_DC      screenDCPin
#define TFT_CS      screenCSPin
#define TFT_RST     255  // 255 = unused, connect to 3.3V
#define TFT_MOSI    sdiPin
#define TFT_SCLK    sckPin
#define TFT_MISO    sdoPin
ILI9341_t3 tft = ILI9341_t3(TFT_CS, TFT_DC, TFT_RST, TFT_MOSI, TFT_SCLK, TFT_MISO);

#define width 320
#define height 240
#define menuWidth 140
#define menuItemHeight 48
#define statusWidth (width-menuWidth)
#define statusHeight 32
#define navButtonWidth (statusWidth/2)
#define visualizerWidth (statusWidth)
#define visualizerHeight (height - (menuItemHeight * 2 + statusHeight * 2))

#define screenMenuLen 12

struct Point {
  Point() {
    x = 0;
    y=0;
  }
  Point(uint16_t x, uint16_t y): x{x}, y{y} {
  };
  uint16_t x;
  uint16_t y;
};

struct Rectangle {
  Rectangle() {
    p1 = Point();
    p2 = Point();
  }
  Rectangle(Point p1, Point p2) : p1{p1}, p2{p2} {};
  Rectangle(uint16_t x1, uint16_t y1, uint16_t x2, uint16_t y2) {
    p1 = Point(x1, y1);
    p2 = Point(x2, y2);
  }

  Point p1;
  Point p2;
};

struct Window {
  Rectangle extent;
  String text;
  uint16_t bgcolor;
  uint16_t fgcolor;
  bool enabled;
  bool redraw;
  bool highlight;
  uint32_t textUsecs;
};

#define numWindows 12

enum windowIndex {
  menuText1 = 0,
  menuText2,
  menuText3,
  menuText4,
  menuText5,
  backText,
  fwdText,
  cancelText,
  okText,
  statusBar1,
  visualizerWindow,
  statusBar2
};


struct Window windows[numWindows];

struct Point cursor(const struct Window &window) {
  uint16_t left = window.extent.p1.x;
  //uint16_t right = window.extent.p2.x;
  uint16_t top = window.extent.p1.y;
  uint16_t bottom = window.extent.p2.y;
  uint16_t bheight = bottom-top;

  return Point(left+5, top+(bheight/2) - 9);
}

void setWindowCursor(const struct Window &window) {
  Point c = cursor(window);
  tft.setCursor(c.x, c.y);
}

uint32_t brightness = 127;
uint32_t brightnessSet = 255;

bool moreUp = false;
bool moreDown = false;

void renderScreen(uint32_t deltaUsecs) {
  for (int i=0; i<numWindows; i++) {
    if (!windows[i].redraw) {
      continue;
    }

    Rectangle r = windows[i].extent;
    if (!windows[i].enabled) {
      tft.fillRect(r.p1.x+1, r.p1.y+1, (r.p2.x-r.p1.x)-2, (r.p2.y-r.p1.y)-2, ILI9341_BLACK);
      windows[i].redraw = false;
      continue;
    }

    uint16_t fgcolor, bgcolor;

    if (windows[i].highlight) {
      fgcolor = windows[i].bgcolor;
      bgcolor = windows[i].fgcolor;
    } else {
      fgcolor = windows[i].fgcolor;
      bgcolor = windows[i].bgcolor;
    }

    tft.fillRect(r.p1.x+1, r.p1.y+1, (r.p2.x-r.p1.x)-2, (r.p2.y-r.p1.y)-2, bgcolor);
    
    tft.setTextColor(fgcolor);
    tft.setTextSize(2);
    setWindowCursor(windows[i]);
    tft.println(windows[i].text);

    windows[i].redraw = false;
  }

  if (moreUp) {
    tft.setTextColor(windows[menuText1].highlight ? windows[menuText1].bgcolor : windows[menuText1].fgcolor);
    tft.setCursor(52, 0);
    tft.println("...");
  }

  if (moreDown) {
    tft.setTextColor(windows[menuText5].highlight ? windows[menuText5].bgcolor : windows[menuText5].fgcolor);
    tft.setCursor(52, 218);
    tft.println("...");
  }

  if (brightness != brightnessSet) {
    analogWrite(backlightPin, brightness);
    brightnessSet = brightness;
    Serial.println("set brightness to " + String(brightness) + "/255");
  }
}

void screenSetup() {
  Serial.println(" 1");
  for (int i=0; i<numWindows; i++) {
    windows[i].extent = Rectangle(0,0,0,0);
    windows[i].text = "";
    windows[i].bgcolor = ILI9341_DARKGREY;
    windows[i].fgcolor = ILI9341_WHITE;
    windows[i].enabled = true;
    windows[i].redraw = true;
    windows[i].highlight = false;
    windows[i].textUsecs = 0;
  }


  for (int i=0; i<5; i++) {
    windows[i].extent = Rectangle (0, i*menuItemHeight, menuWidth, (i+1)*menuItemHeight);
  }

  Serial.println("2");
  windows[backText].extent  = Rectangle (menuWidth, 0, menuWidth + navButtonWidth, menuItemHeight);
  windows[fwdText].extent  = Rectangle (menuWidth + navButtonWidth, 0, width, menuItemHeight);
  windows[cancelText].extent  = Rectangle (menuWidth, menuItemHeight, menuWidth + navButtonWidth, menuItemHeight * 2);
  windows[okText].extent  = Rectangle (menuWidth + navButtonWidth, menuItemHeight, width, menuItemHeight * 2);

  windows[statusBar1].extent = Rectangle (menuWidth, menuItemHeight * 2, width, menuItemHeight * 2 + statusHeight);
  windows[statusBar1].bgcolor = ILI9341_BLACK;
  windows[statusBar2].extent  = Rectangle (menuWidth, height-statusHeight, width, height);
  windows[statusBar2].bgcolor = ILI9341_BLACK;

  windows[visualizerWindow].extent = Rectangle (menuWidth, menuItemHeight * 2 + statusHeight, width, height - statusHeight);
  windows[visualizerWindow].bgcolor = ILI9341_BLACK;

  windows[backText].text = "back";
  windows[backText].enabled = false;
  windows[fwdText].text = "forward";
  windows[fwdText].enabled = false;
  windows[cancelText].text = "cancel";
  windows[cancelText].enabled = false;
  windows[okText].text = "ok";
  windows[okText].enabled = false;

  Serial.println("3");
  pinMode(backlightPin, OUTPUT);
  analogWriteFrequency(backlightPin, 3611*2); /* default is 3.611 kHz */
  analogWrite(backlightPin, brightness);
  brightnessSet = brightness;

  Serial.println("3.1 tft:" + String((uint32_t)&tft));
  delayMicroseconds(100000);

  tft.begin();
  Serial.println("3.2");
  delayMicroseconds(100000);
  tft.setRotation(1);
  tft.setClock(100000000);
  tft.fillScreen(ILI9341_BLACK);

  Serial.println("4");
  renderScreen(0);

  uint8_t x = tft.readcommand8(ILI9341_RDMODE);
  Serial.print("Display Power Mode: 0x"); Serial.println(x, HEX);
  x = tft.readcommand8(ILI9341_RDMADCTL);
  Serial.print("MADCTL Mode: 0x"); Serial.println(x, HEX);
  x = tft.readcommand8(ILI9341_RDPIXFMT);
  Serial.print("Pixel Format: 0x"); Serial.println(x, HEX);
  x = tft.readcommand8(ILI9341_RDIMGFMT);
  Serial.print("Image Format: 0x"); Serial.println(x, HEX);
  x = tft.readcommand8(ILI9341_RDSELFDIAG);
  Serial.print("Self Diagnostic: 0x"); Serial.println(x, HEX); 
  Serial.println("5");
}

/*
 * Update one column of visualizer -- typically the height of each column will
 * correspond to the value from that analog input.
 * We can display 3 values at once by treating red green and blue separately.
 * r, g, and b should be from 0.0 to 1.0.
 */
void visualizerUpdateGraph(int column, float r, float g, float b) {
  if (windows[visualizerWindow].enabled == true) {
    return;
  }

  if (column < 0 || column >= visualizerWidth) {
    return;
  }

  r = clamp(r);
  g = clamp(g);
  b = clamp(b);

  int rHeight = r * (visualizerHeight - 4);
  int gHeight = g * (visualizerHeight - 4);
  int bHeight = b * (visualizerHeight - 4);

  int minHeight = 0;
  int midHeight = 0;
  int maxHeight = 0;

  uint16_t cmin = ILI9341_WHITE;
  uint16_t c2 = 0;
  uint16_t c3 = 0;
  uint16_t cmax = ILI9341_BLACK;

  if (r <= g && r <= b) {
    minHeight = rHeight;
    c2 = ILI9341_CYAN;
    if (g <= b) {
      midHeight = gHeight;
      maxHeight = bHeight;
      c3 = ILI9341_BLUE;
    } else {
      midHeight = bHeight;
      maxHeight = gHeight;
      c3 = ILI9341_GREEN;
    }
  } else if (g <= r && g <= b) {
    minHeight = gHeight;
    c2 = ILI9341_MAGENTA;
    if (r <= b) {
      midHeight = rHeight;
      maxHeight = bHeight;
      c3 = ILI9341_BLUE;
    } else {
      midHeight = bHeight;
      maxHeight = rHeight;
      c3 = ILI9341_RED;
    }
  } else {
    minHeight = bHeight;
    c2 = ILI9341_YELLOW;
    if (r <= g) {
      midHeight = rHeight;
      maxHeight = gHeight;
      c3 = ILI9341_GREEN; 
    } else {
      midHeight = gHeight;
      maxHeight = rHeight;
      c3 = ILI9341_RED;
    }
  }

  int x = windows[visualizerWindow].extent.p1.x + column + 1;
  int y = windows[visualizerWindow].extent.p2.y - 2;
  tft.drawLine(x, y, x, y-minHeight, cmin);
  tft.drawLine(x, y-(minHeight+1), x, y-midHeight, c2);
  tft.drawLine(x, y-(midHeight+1), x, y-maxHeight, c3);
  tft.drawLine(x, y-(maxHeight+1), x, y-(visualizerHeight-4), cmax);
}

// MENU

enum menuItemType {
  action,
  toggle,
  value,
  selection,
  submenu,
  empty
};


/* menu item that's selected for editing */
struct MenuItem* editItem = nullptr;

struct MenuItem {
  MenuItem(String text, menuItemType type,
          struct MenuItem* m1 = nullptr,
          struct MenuItem *m2 = nullptr,
          struct MenuItem *m3 = nullptr,
          struct MenuItem *m4 = nullptr,
          struct MenuItem *m5 = nullptr) : text{text}, type{type} {
    select = nullptr;
    numChildren = 0;
    data = 0;

    if (m1 != nullptr && type == submenu) {
      children[numChildren++] = m1;
      if (m2 != nullptr) {
        children[numChildren++] = m2;
        if (m3 != nullptr) {
          children[numChildren++] = m3;
          if (m4 != nullptr) {
            children[numChildren++] = m4;
            if (m5 != nullptr) {
              children[numChildren++] = m5;  
            }
          }
        }
      }
    }
  }
  MenuItem(String text, menuItemType type, struct MenuItem **items, uint16_t numItems) : text{text}, type{type}, childrenExtended{items}, numChildren{numItems} {
    for (int i = 0; i < 5 && i < numItems; i++) {
      children[i] = items[i];
    }
  }
  MenuItem(String text, void (*select)(void *data), void *data = nullptr) : text{text}, type{action}, select{select}, data{data} {}
  MenuItem(String text, enum menuItemType type, void *data, uint32_t *minValue, uint32_t *maxValue) : text{text}, type{type}, data{data}, minValue{minValue}, maxValue{maxValue} {}
  MenuItem(String text, enum menuItemType type, void *data, uint32_t defaultValue = 0, void (*select)(void *data) = nullptr) : text{text}, type{type}, select{select}, data{data}, defaultValue{defaultValue} {
    if (type == toggle && data != nullptr) {
      highlight = *((bool*)data);
    }

    if (type == selection) {
      highlight = defaultValue == *(uint32_t*)data;
    }
  }


  bool checkHighlight() {
    switch (type) {
      case toggle:
        highlight = *(bool*)data != false;
        break;
      case action:
        highlight = false;
        break;
      case selection:
        highlight = *(uint32_t*)data == defaultValue;
        break;
      case value:
        highlight = editItem == this;
        break;
      default:
        highlight = false;
        break;
    }

    return highlight;
  }

  String text;
  enum menuItemType type;
  void (*select)(void* data) = nullptr;
  void *data = nullptr;
  uint32_t defaultValue = 0; /* for "selection" type, the value to set data to */
  uint32_t *minValue = nullptr;
  uint32_t *maxValue = nullptr;
  struct MenuItem** childrenExtended = nullptr;
  uint16_t numChildren = 0;
  struct MenuItem* children[5];
  int scrollOffset = 0;
  bool highlight = false;
};

#define menuStackSize 10
struct MenuItem* menuStack[menuStackSize];
uint16_t menuStackPos = 0; /* points to first unoccupied slot */

struct MenuItem emptyMenuItem = MenuItem("", empty);
struct MenuItem* menu[5] = {&emptyMenuItem, &emptyMenuItem, &emptyMenuItem, &emptyMenuItem, &emptyMenuItem};

void statusTextUpdate();

void menuSelect(struct MenuItem *item, uint16_t button) {
  if (item->type == submenu) {
    for (int i=0; i < 5; i++) {
      if (i + item->scrollOffset < item->numChildren) {
        windows[i].bgcolor = ILI9341_DARKGREY;
        windows[i].text = item->children[i]->text;
        windows[i].enabled = item->children[i]->type != empty;
        menu[i] = item->children[i];
        Serial.println("displaying menu item " + String(i) + " " + windows[i].text);
      } else {
        windows[i].enabled = false;
        menu[i] = &emptyMenuItem;
      }
      windows[i].redraw = true;
    }

    windows[backText].enabled = menuStackPos > 0;
    windows[backText].redraw = true;

    menuStack[menuStackPos] = item;
    menuStackPos++;

    if (menuStackPos >= menuStackSize) {
      menuStackPos = menuStackSize - 1;
      Serial.println("menu stack overflow");
    }

    moreUp = item->scrollOffset > 0;
    moreDown = item->scrollOffset + 5 < item->numChildren; 
  }

  if (item->type == toggle){
    if (item->data != nullptr) {
      *((bool*)(item->data)) = !*(bool*)(item->data);
    }
  }

  if (item->type == selection) {
    if (item->data != nullptr) {
      *((uint32_t*)(item->data)) = item->defaultValue;
    }
  }

  if (item->type == value) {
    windows[visualizerWindow].text = String(*(uint32_t *)(item->data));
    windows[visualizerWindow].enabled = true;
    editItem = item;
    windows[visualizerWindow].redraw = true;
  } else if (windows[visualizerWindow].enabled == true) {
    editItem = nullptr;
    windows[visualizerWindow].text = "";
    windows[visualizerWindow].redraw = true;
    windows[visualizerWindow].enabled = false;
  }

  if (item->select != nullptr) {
    item->select(item->data);
  }

  if (menuStackPos > 0) {
    struct MenuItem *parent = menuStack[menuStackPos-1];
    for (int i = 0; i + parent->scrollOffset < parent->numChildren && i < 5; i++) {
      auto child = parent->children[i];
      bool prevHighlight = child->highlight;

      windows[i].highlight = child->checkHighlight();
      if (windows[i].highlight != prevHighlight) {
        windows[i].redraw = true;
      }
    }
  }

  if (item->type != submenu) {
    statusTextUpdate();
    windows[statusBar1].redraw = true;
    windows[statusBar2].redraw = true;
  }
}

void menuScroll(int offset) {
  struct MenuItem *item = menuStack[menuStackPos-1];

  if (item == nullptr || item->type != submenu || item->childrenExtended == nullptr) {
    return;
  }

  if (offset > 0) {
    if (item->scrollOffset + offset + 5 > item->numChildren) {
      offset = 0;
    }
  } else {
    if (item->scrollOffset + offset < 0) {
      offset = 0;
    }
  }

  item->scrollOffset += offset;

  Serial.println("menuScroll offset " + String(offset));

  moreUp = item->scrollOffset > 0;
  moreDown = item->scrollOffset + 5 < item->numChildren;

  for (int i = 0; i < 5; i++) {
    if (i + item->scrollOffset < item->numChildren) {
      item->children[i] = item->childrenExtended[i + item->scrollOffset];

      Serial.println("menu " + String(i) + " " + item->children[i]->text);
    } else {
      item->children[i] = nullptr;
    }

    if (offset != 0) {
      menu[i] = item->children[i] == nullptr ? &emptyMenuItem : item->children[i];
      windows[i].text = menu[i]->text;
      windows[i].highlight = item->children[i]->checkHighlight();
      windows[i].redraw = true;
    }
  }
}

void menuPress(uint8_t button, uint16_t pressure, uint32_t deltaUsecs) {
  windows[button].highlight = !windows[button].highlight;
  windows[button].redraw = true;
}

void menuRelease(uint8_t button, uint16_t pressure, uint32_t deltaUsecs) {
  windows[button].highlight = !windows[button].highlight;
  windows[button].redraw = true;
  if (button <= menuText5) {
    Serial.println("released button " + String(button));
    menuSelect(menu[button], button);
  } else if (button == backText) {
    Serial.println("back button released");
    menuBack();
  }
}

void menuBack() {
  if (menuStackPos > 1) {
    menuStackPos-=2;

    windows[backText].enabled = menuStackPos > 0;
    windows[backText].redraw = true;
    menuSelect(menuStack[menuStackPos], menuStackPos);
  }
}

void menuUpdate(uint32_t deltaUsecs) {
  struct Window* status1 = &windows[statusBar1];
  if (status1->textUsecs > 0) {
    if (deltaUsecs >= status1->textUsecs) {
      status1->textUsecs = 0;
      status1->text = "";
      status1->redraw = true;
    } else {
      status1->textUsecs -= deltaUsecs;
    }
  }
}

/* CAN BUS */

FlexCAN_T4<CAN1, RX_SIZE_256, TX_SIZE_16> can;

void canSetup() {
  can.begin();
  can.setBaudRate(1000000);
}

void canUpdate() {
  CAN_message_t msg;
  if (can.read(msg)){
    Serial.print("can message received ID 0x");
    Serial.print(msg.id, HEX);
    Serial.print(" data ");
    for (int i = 0; i < 8; i++) {
      Serial.print(msg.buf[i], HEX); Serial.print(" ");
    }
    Serial.println();
  }

  /*
  static uint32_t t_start = millis();
  if (millis() > t_start + 100) {
    Serial.println("sending ping");
    msg.id = random(0x1,0x7FE);
    for ( uint8_t i = 0; i < 8; i++ ) msg.buf[i] = i + 1;
    can.write(msg);
    t_start = millis();
  } */
}

/* MIDI */

bool useUsbMidi = true;
bool useDinMidi = false;

int noteOnCount = 0;
int noteOffCount = 0;

struct MidiLocalSettings : public MIDI_NAMESPACE::DefaultSettings {
  static const bool UseRunningStatus = true;  /* avoid re-sending status byte when it hasn't changed */
  static const unsigned SysExMaxSize = 10;     /* we don't expect to receive sysex */
};

MIDI_CREATE_CUSTOM_INSTANCE(HardwareSerial, Serial5, dinMidi, MidiLocalSettings);

//MIDI_CREATE_INSTANCE(HardwareSerial, Serial5, dinMidi);

uint64_t midiMsgsSent = 0;
uint64_t midiMsgsReceived = 0;
int pressureBackoff = 5000;

#define doMidi(func, ...) { \
  midiMsgsSent++; \
  if (useUsbMidi) { \
    usbMIDI.func(__VA_ARGS__); \
  } \
  if (useDinMidi) { \
    dinMidi.func(__VA_ARGS__); \
  } \
}

#define doMidiDivergent(usbFunc, serialFunc, ...) { \
  midiMsgsSent++; \
  if (useUsbMidi) { \
    usbMIDI.usbFunc(__VA_ARGS__); \
  } \
  if (useDinMidi) { \
    dinMidi.serialFunc(__VA_ARGS__); \
  } \
}

void midiNoteOn(uint8_t note, uint8_t velocity, uint8_t channel) {
  doMidi(sendNoteOn, note, velocity, channel);
  noteOnCount++;
}

void midiNoteOff(uint8_t note, uint8_t velocity, uint8_t channel) {
  doMidi(sendNoteOff, note, velocity, channel);
  noteOffCount++;
}

void midiNoteOffNoCount(uint8_t note, uint8_t velocity, uint8_t channel) {
  doMidi(sendNoteOff, note, velocity, channel);
}

void midiPitchBend(int16_t pb, uint8_t channel) {
  doMidi(sendPitchBend, pb, channel);
}

void midiAfterTouch(uint8_t volume, uint8_t channel) {
  doMidi(sendAfterTouch, volume, channel);
}

void midiPolyAfterTouch(uint8_t note, uint8_t pressure, uint8_t channel) {
  // sendPolyPressure provokes a deprecation warning,
  // but sendAfterTouch in the usb midi api doesn't have
  // a polyphonic variant
  doMidiDivergent(sendAfterTouchPoly, sendAfterTouch, note, pressure, channel);
}

void midiControlChange(uint8_t cc, uint8_t value, uint8_t channel) {
  if (cc > 127) {
    return;
  }
  doMidi(sendControlChange, cc, value, channel);
}

void midiProgramChange(uint8_t bank, uint8_t channel) {
  doMidi(sendProgramChange, bank, channel);
}

void midiRPN14Bit(uint8_t pmsb, uint8_t plsb, uint8_t vmsb, uint8_t vlsb, uint8_t channel) {
  midiControlChange(0x65, pmsb, channel);
  midiControlChange(0x64, plsb, channel);
  midiControlChange(0x6,  vmsb, channel);
  midiControlChange(0x26, vlsb, channel);
  midiControlChange(0x65, 127, channel);
  midiControlChange(0x64, 127, channel);
}

void midiRPN(uint8_t pmsb, uint8_t plsb, uint8_t v, uint8_t channel) {
  midiControlChange(0x65, pmsb, channel);
  midiControlChange(0x64, plsb, channel);
  midiControlChange(0x6,  v, channel);
  midiControlChange(0x65, 127, channel);
  midiControlChange(0x64, 127, channel);
}

void midiSysEx(uint32_t len, const byte *in) {
  doMidi(sendSysEx, len, in, true);
}

void sendDummySysEx() {
  int len = 100;

  byte buf[len];

  for (int i=0; i<len; i++) {
    buf[i] = i % 11;
  }

  midiSysEx(len, buf);
  Serial.println("sent dummy sysex");
}

int midiBufferSize = 0;

void midiSetup(){
  Serial.println("serial fifo size " + String(Serial5.availableForWrite()));
  dinMidi.begin();
  midiBufferSize = Serial5.availableForWrite();
  Serial.println("midiBufferSize set to " + String(midiBufferSize)); /* the default size appears to be 39 bytes */
}

int midiBufferInUse() {
  if (useDinMidi) {
    int avail = Serial5.availableForWrite();
    int inUse = midiBufferSize - avail;
    //Serial.println("midiBufferInUse " + String(inUse) + "/" + String(midiBufferSize) + " : " + String(avail));
    if (inUse < 0) {
      return 0;
    }
    return inUse;
  }
  return 0;
}

bool midiReady() {
  int inUse = midiBufferInUse();
  if (inUse > 1000 || inUse < -1000) {
    Serial.println("midiReady inUse = " + String(inUse));
  }

  return midiBufferInUse() < 15;
}

bool midiReadyLowPriority() {
  return midiBufferInUse() < 10;
}

bool midiIdle() {
  return midiBufferInUse() == 0;
}

void midiReadyWait() {
  int iterations = 0;
  while (!midiReady()) {
    
    delayMicroseconds(5);
    iterations++;
    if (iterations+1 % 10 == 0) {
      Serial.println("midiReadyWait " + String(iterations) + " " + String(midiBufferInUse()));
    }
  };
}


float mpePolyAfterTouchMin = 0.75f;
float mpePolyAfterTouchMax = 1.25f;


/* MPE */

int mpeChannels = 16;
int mpeChannelsSet = mpeChannels;
int firstMpeChannel = 0;

#define noOne 0xffff

struct MpeChannelState{
  int channel;
  bool playing;
  uint32_t age;
  uint8_t lastNote;
  uint8_t lastVolume;
  uint8_t lastPolyAT;
  float volume;
  uint8_t lastFilter;
  int16_t lastPitchBend;
  double lastBendInterval;
  uint32_t pitchBendAge;
  uint16_t owner;
  uint8_t lastProgramChangeSent;
  uint32_t volumeAge;
  uint8_t lastBankMsb;
  uint8_t lastBankLsb;
  double originalPitch;
  void (*stealCallback)(uint16_t owner);
};

/* Store state of "global" CCs that should affect all channels the same. */
struct CCState {
  CCState() {};
  uint16_t dirty = 0;       /* bitmask of channels requiring CC updates */
  uint8_t next = 127;       /* keep track of next in-use CC so we can iterate faster */
  uint8_t value = 128;      /* curent value of CC, or 128 if either unused or not set */
  bool noMulticast = 0;     /* do not send on all channels, as we control those individually elsewhere */
  bool forceMulticast = 0;  /* multicast even on MPE */
  uint8_t center = 0;       /* "neutral" value */
};

/* most of these will go unused, but that's okay */
struct CCState ccs[128];

struct MpeChannelState mpeState[16];
uint32_t programChange = 0;

bool doProgramChange = false;
bool doMsb = false;
bool doLsb = false;

void mpeMulticastCC(uint8_t cc, float value) {
  int v = (int)(value * 127.0f);

  Serial.println("mpeMulticastCC " + String(cc) + " " + String(v));

  if (v > 127) {
    v = 127;
  } else if (v < 0) {
    v = 0;
  }

  if (v == ccs[cc].value) {
    return;
  }

  ccs[cc].value = v;
  ccs[cc].dirty = 0xffff;

  for (int i = cc-1; i >= 0; i--) {
    if (ccs[i].next <= cc) {
      break;
    }

    ccs[i].next = cc;
  }
}

uint32_t mpeIdleScore(int channel) {
  struct MpeChannelState *state = &mpeState[channel];
  uint32_t score = 0;

  if (state->playing == false) {
    score = 0x80000000;
  }

  if (state->age >= 0x80000000) {
    score += 0x7fffffff;
  } else {
    score += state->age;
  }
  return score;
}

bool skipChannel10 = false;

struct MpeChannelState *getMpeChannel() {
  int bestChannel = 0;
  uint32_t bestScore = 0;
  if (!midiReady()) {
    return nullptr;
  }
  for (int channel = firstMpeChannel; channel < firstMpeChannel+mpeChannels; channel++) {
    if (channel+1 == 10 && skipChannel10) {
      continue;
    }

    uint32_t score = mpeIdleScore(channel);
    if (score > bestScore) {
      bestScore = score;
      bestChannel = channel;
    }
  }

  struct MpeChannelState *state = &mpeState[bestChannel];
  if (state->playing  && state->age > 100000) {
    /* steal note */
    Serial.println("stealing note");
    midiNoteOff(state->lastNote, 127, bestChannel+1);
    state->playing = false;
    if (state->stealCallback != nullptr) {
      state->stealCallback(state->owner);
    }
    return state;
  }

  if (state->playing) {
    return nullptr;
  }

  return state;
}

#define middleC 60 // midi note number for C4

uint32_t pbRange = 2;  // could use float, but int is more compatible with MIDI pitch bend range configuration
uint32_t pressureCC = 128;
uint32_t pressureCCSet = pressureCC;
uint32_t maxMpePressure = 127;

/* seconds from full to none, or vice versa */
double attack = 0.0;
double decay = 1.0;

double pitchToCents(double pitch) {
  return (log(pitch) / log(2.0)) * 1200.0;
}

double pitchToOctaves(double pitch) {
  return (log(pitch) / log(2.0));
}

double centsToPitch(double cents) {
  const double factor = pow(2.0, (1.0/1200.0));
  return pow(factor, cents);
}

double semitone = pow(2.0, (1.0/12.0));

enum midiType {
  monovoice,
  monotimbral,
  tuningtable,
  multitimbral,
  mpe
};

enum tuningTableType {
  noTuningTable,
  eTuningTable,
  mtsTuningTable
};

enum tuningTableType tuningTableType;

bool doMpeDynamicVelocity = true;
bool doMpeDynamicPressure = true;
bool doMpePolyAfterTouch = false;
float mpeStaticVelocity = 0.75;
float minVelocity = 0.2;

double transpose = 1.0;
uint32_t mpeBankLsbMin = 0;
uint32_t mpeBankLsbMax = 0;
uint32_t mpeBankMsbMin = 0;
uint32_t mpeBankMsbMax = 0;
uint32_t mpeBankMsb = 0;
uint32_t mpeBankLsb = 0;
bool unlockBankRange = false;
bool unlockBankRangeSet = unlockBankRange;
enum midiType midiType = mpe;

double masterPbUpRange = 3.0/2.0;
double masterPbDownRange = 2.0/3.0;

double pbUp = 0.0;
double pbDown = 0.0;

uint32_t masterPbAge = 0;
int16_t lastMasterPb = 0;
uint32_t masterPbBackoff = 8000;

double masterPbRange = 12.0;

bool forcePbRange = false;

// 1/seconds to fall from full intensity
float releaseRate = 2.0f;

float pressureExponent = 0.8f;

bool delayNoteOff = false;
bool noteOnFirst = false;

double pitchReferenceHz() {
  double c = 440.0 / (pow(pow(2, 1.0/12.0), 9));

  return c * transpose;
}

int16_t calculatePitchBend(double pbUp, double pbDown, double pitch, double pbRange) {
  double pb = clampDoubleBipolar(pbUp - pbDown);

  double cents = pb > 0
    ? pitchToCents(masterPbUpRange) * pb
    : pitchToCents(masterPbDownRange) * -pb;

  cents += pitchToCents(pitch);
    
  double range = (double)MIDI_PITCHBEND_MAX - (double)MIDI_PITCHBEND_MIN;
  int pbi = (int)((range / (pbRange * 2.0)) * (cents / 100.0));
  if (pbi > MIDI_PITCHBEND_MAX) {
    return MIDI_PITCHBEND_MAX;
  } else if (pbi < MIDI_PITCHBEND_MIN) {
    return MIDI_PITCHBEND_MIN;
  }

  return pbi;
}

void prepareChannel(struct MpeChannelState *state) {
  int channel = state->channel;

  midiReadyWait();

  /*
   * If we're updating the LSB, we have to send the MSB first
   * even if it hasn't changed. (On the XV-2020 at least.)
   */
  if (doMsb && (state->lastBankMsb != mpeBankMsb || state->lastBankLsb != mpeBankLsb)) {
    midiControlChange(0, mpeBankMsb, channel + 1);
    state->lastBankMsb = mpeBankMsb;
    state->lastProgramChangeSent = 128;  /* force new program change to be sent */
    Serial.println("sent msb " + String(mpeBankMsb) + " channel " + String(channel));
  }

  if (doLsb && state->lastBankLsb != mpeBankLsb) {
    midiControlChange(32, mpeBankLsb, channel + 1);
    state->lastBankLsb = mpeBankLsb;
    state->lastProgramChangeSent = 128; /* force new program change to be sent */
    Serial.println("sent lsb" + String(mpeBankLsb) + " channel " + String(channel));
  }

  if (doProgramChange && state->lastProgramChangeSent != programChange) {
    midiProgramChange(programChange, channel + 1);
    state->lastProgramChangeSent = programChange;

    // set pitch bend range
    if (forcePbRange) {
      midiReadyWait();
      midiRPN(0, 0, pbRange, state->channel+1);
    }
    Serial.println("sent program change " + String(programChange) + " channel " + String(channel));
  }
}

bool bendUpOnly = false;

struct MpeChannelState *beginMpeNote(double pitch, double velocity, double pressure, uint16_t owner, void stealCallback(uint16_t)) {

  int v = (int) (1.0 + (velocity*126));
  if (v > 127) {
    v = 127;
  } else if (v < 1) {
    v = 1;
  }

  double shift = pitchToCents(pitch * transpose);
  //state->originalPitch = shift;

  int semitones = shift / 100.0f;
  float cents = shift - (semitones * 100.0f);
  if (cents > 50.0f) {
    semitones++;
    cents -= 100.0f;
  } else if (cents < -50.0f) {
    semitones--;
    cents += 100.0f;
  }

  if (bendUpOnly && cents < 0.0f) {
    semitones--;
    cents += 100.0f;

  }

  int note = middleC + semitones;
  int16_t pb;
  if (midiType == mpe) {
    pb = calculatePitchBend(0.0, 0.0, centsToPitch(cents), (double)pbRange);
  } else {
    pb = calculatePitchBend(pbUp, pbDown, centsToPitch(cents), (double)pbRange);
  }

  if (note < 0) {
    Serial.print("out of range note (too low) ");
    Serial.println(note);
    return nullptr;
  } else if (note > 127) {
    Serial.print("out of range note (too high) ");
    Serial.println(note);
    return nullptr;
  }


  struct MpeChannelState *state = getMpeChannel();
  if (state == nullptr) {
    Serial.println("no available channels");
    return nullptr;
  }

  if (state->playing) {
    Serial.println("getMpeChannel returned active channel (shouldn't happen)");
    return nullptr;
  }

  state->age = 0;
  state->owner = owner;
  state->stealCallback = stealCallback;

  if(midiReady()) {
    state->playing = true;
    state->lastNote = note;
    state->lastBendInterval = centsToPitch(cents);

    if (midiType != mpe) {
      prepareChannel(state);
    }

    midiPitchBend(pb, state->channel+1);
    state->lastPitchBend = pb;
    state->volumeAge = 1000000;
    state->volume = 0.0f; /* we'll reset this in continueMpeNote, but we need to zero it first */

    if (noteOnFirst) {
      midiNoteOn(note, v, state->channel+1);
    }

    if (!continueMpeNote(state, pressure, 0)) {
      state->playing = false;
      state->owner = noOne;
      return nullptr;
    }

    if (!noteOnFirst) {
      midiNoteOn(note, v, state->channel+1);
    }

    Serial.println("note-on channel" + String(state->channel+1) + " note " + String(note) + " cents " + String(cents) + " (" + String(centsToPitch(cents)) + ") pb " + String(pb) + " pressure " + pressure);
  }

  state->age = 0;

  return state;
}

bool continueMpeNote(struct MpeChannelState *state, double pressure, uint32_t deltaUsecs) {
  uint32_t concurrentNotes = noteOnCount - noteOffCount;

  float minPressure = state->volume - releaseRate * ((float)deltaUsecs/1000000.0f);
  if (pressure < minPressure) {
    pressure = minPressure;
  }

  if (pressure < 0.0f) {
    pressure = 0.0f;
  }

  state->volume = pressure;

  /* rate limit pressure updates */
  if (state->volumeAge + deltaUsecs < pressureBackoff * concurrentNotes) {
    state->volumeAge += deltaUsecs;
    return true;
  }

  int max = maxMpePressure;

  uint16_t volume = doMpeDynamicPressure ? pressure * (float)max : max;
  if (volume > max) {
    volume = max;
  }

  if(midiReadyLowPriority()) {
    if (volume != state->lastVolume) {
      if (midiType == mpe) {
        midiAfterTouch(volume, state->channel+1);
        //Serial.println("aftertouch " + String(volume));
      } else {
        if (pressureCC != 128) {
          midiControlChange(pressureCC, volume, state->channel+1);
        }
      }
      state->lastVolume = volume;
      state->volumeAge = 0;
    } else {
      //Serial.println("aftertouch " + String(volume) + " (not sent, no change)");
    }

    if (doMpePolyAfterTouch) {
      uint8_t atPressure = pressure <= mpePolyAfterTouchMin
        ? 0
        : pressure >= mpePolyAfterTouchMax
          ? 127
          : ((pressure - mpePolyAfterTouchMin) / (mpePolyAfterTouchMax - mpePolyAfterTouchMin)) * 127;
      if (atPressure != state->lastPolyAT) {
        midiPolyAfterTouch(state->lastNote, atPressure, state->channel+1);
        state->lastPolyAT = atPressure;
        //Serial.println("poly aftertouch " + String(atPressure));
      }
    }
  } else {
    Serial.print(".");
    state->volumeAge += deltaUsecs;
  }

  if (pressure <= 0.0f && delayNoteOff) {
    if (endMpeNote(state)) {
      return false;
    }
  }
  
  return true;
}

bool endMpeNote(struct MpeChannelState *state) {
  if (midiReady()) {
    state->age = 0;
    /* if note wasn't already stolen */
    if (state->playing == true) {
      midiNoteOff(state->lastNote, 63, state->channel+1);
      state->lastPolyAT = 0;
    }
    state->playing = false;
    Serial.println("sent note-off channel " + String(state->channel+1) + " note " + String(state->lastNote));

    return true;
  }
  return false;
}

void doMpeMasterPitchbend(double pbUp, double pbDown, uint32_t deltaUsecs) {
  if (deltaUsecs < masterPbBackoff - masterPbAge) {
    masterPbAge += deltaUsecs;
    return;
  }

  if (midiType == mpe || midiType == tuningtable) {
    int16_t pbi = calculatePitchBend(pbUp, pbDown, 1.0, masterPbRange);
    if (pbi != lastMasterPb) {
      midiPitchBend(pbi, 1);
      lastMasterPb = pbi;
    }
  } else if (midiType == multitimbral || midiType == monotimbral || midiType == monovoice) {
    for (int channel = firstMpeChannel; channel < firstMpeChannel+mpeChannels; channel++) {
      if (mpeState[channel].playing) {
        int16_t pbi = calculatePitchBend(pbUp, pbDown, mpeState[channel].lastBendInterval, (double)pbRange);
        if (mpeState[channel].lastPitchBend != pbi) {
          midiPitchBend(pbi, channel+1);
          mpeState[channel].lastPitchBend = pbi;
        }
      }
    }
  }

  masterPbAge = 0;
}

/* TUNING TABLE MIDI */

struct MidiNote {
  MidiNote(uint8_t noteNumber, uint8_t channel) : noteNumber{noteNumber}, channel{channel} {};
  uint8_t noteNumber = 0;
  uint8_t channel = 0;
};

struct MidiTuningTableEntry {
  uint8_t noteNumber;
  uint8_t channel;
  float pitch;
  uint16_t users;
  uint8_t lastPressure;
  uint32_t pressureAge;
};



struct MidiTuningTableEntry *midiTuningTable = nullptr;
int midiTuningTableSize = 0;

struct MidiTuningTableEntry *tuningTableLookup(float pitch, struct MidiTuningTableEntry tt[], int ttSize) {

  /* tuning table isn't guaranteed to be sorted; we could sort it and speed this up with binary search */
  float bestError = 1000.0f;
  struct MidiTuningTableEntry *bestMatch = nullptr;

  for (int i = 0; i < ttSize; i++) {
    float error = abs((pitch / tt[i].pitch) - 1.0f);
    if (error < bestError) {
      bestError = error;
      bestMatch = &tt[i];
    }
  }

  return bestMatch;
}


struct MidiTuningTableEntry* beginTuningTableNote(double pitch, double velocity, double pressure, uint16_t owner, void stealCallback(uint16_t)) {
  if (midiTuningTable == nullptr || midiTuningTableSize == 0) {
    return nullptr;
  }

  struct MidiTuningTableEntry *tte = tuningTableLookup(pitch, midiTuningTable, midiTuningTableSize);
  
  if (midiReady()) {
    int v = (int) (1.0 + (velocity*126));
    if (v > 127) {
      v = 127;
    } else if (v < 1) {
      v = 1;
    }

    midiNoteOn(tte->noteNumber, v, tte->channel + 1);
    tte->users++;
    tte->lastPressure = 127;
    tte->pressureAge = 0;
    return tte;
  }

  return nullptr;
}


bool continueTuningTableNote(struct MidiTuningTableEntry *tte, double pressure, uint32_t deltaUsecs) {
  uint32_t concurrentNotes = noteOnCount - noteOffCount;

  if (pressure < 0.0f) {
    pressure = 0.0f;
  }

  /* rate limit pressure updates */
  if (tte->pressureAge + deltaUsecs < pressureBackoff * concurrentNotes) {
    tte->pressureAge += deltaUsecs;
    return true;
  }

  if (!midiReady()) {
    return true;
  }

  if (doMpePolyAfterTouch) {
    int max = mpePolyAfterTouchMax;
    int min = mpePolyAfterTouchMin;
    uint8_t atPressure = pressure <= min
      ? 0
      : pressure >= max
        ? 127
        : ((pressure - min) / (max - min)) * 127;

    if (atPressure != tte->lastPressure) {
      midiPolyAfterTouch(tte->noteNumber, atPressure, tte->channel+1);
      tte->lastPressure = atPressure;
    }

    tte->pressureAge = 0;
  }

  return true; 
}

bool endTuningTableNote(struct MidiTuningTableEntry* tte) {
  if (midiReady()) {
    midiNoteOff(tte->noteNumber, 63, tte->channel + 1);
    tte->users--;
    return true;
  }

  return false;
}

/* return the buffer offset for a given midi note and sub-byte */
#define eTuningTableOffset(index, field) (7 + (index * 4) + field)

/*
 * tuning table setup for the Grey Matter E! expansion board for Yamaha DX7
 */
void sendETuningTable(uint8_t channel = 0) {
  uint8_t buffer[521];

  buffer[0] = 0xf0; // SysEx begin
  buffer[1] = 0x12; // Grey Matter ID number
  buffer[2] = channel & 0xf;
  buffer[3] = 0x00; // 00 = E! for DX7
  buffer[4] = 0x00; // module number (if more than one E! is listening on same midi channel)
  buffer[5] = 0x05; // bank
  buffer[6] = 0x02; // tuning table (1 scale)

  for (int midiKey = 0; midiKey < 128; midiKey++) {
    for (int field = 0; field < 4; field++) {
      buffer[eTuningTableOffset(midiKey, field)] = 0x00;
    }
  }

  for (int i = 0; i < midiTuningTableSize; i++) {
    struct MidiTuningTableEntry *tte = &midiTuningTable[i];
    uint8_t note = tte->noteNumber;

    float pitch = tte->pitch * transpose;
    float octaves = pitchToOctaves(pitch) + 5.0f;
    uint8_t octave = (uint8_t)octaves;

    uint32_t remainder = (octaves - octave) * 4096.0f;

    buffer[eTuningTableOffset(note, 0)] = octave;
    buffer[eTuningTableOffset(note, 1)] = (remainder & 0x0f00) >> 8;
    buffer[eTuningTableOffset(note, 2)] = (remainder & 0x00f0) >> 4;
    buffer[eTuningTableOffset(note, 3)] =  remainder & 0x000f;

    Serial.println("Tuning Table Entry " + String(note) + " octave " + String(octave) + " + " + String(remainder) + " / 4096");
  }

  /* incorrect checksum generates a warning */
  uint8_t sum = 0;
  for (int i = 7; i < 519; i++) {
    sum += buffer[i];
  }

  uint8_t checksum = (((~sum) & 0x7f) + 1) & 0x7f;

  buffer[519] = checksum;
  buffer[520] = 0xf7; /* end SysEx */

  midiSysEx(521, buffer);
}

/* GENERIC VOICE API */

struct VoiceHandle {
  int stolen(int owner) {
    if (this->midiType == tuningtable) {
      return false;
    }

    if (!mpeState) {
      return true;
    }

    return mpeState->owner != owner;
  }

  enum midiType midiType;
  union {
    struct MpeChannelState *mpeState;
    struct MidiTuningTableEntry *tte;
  };
};

bool beginNote(struct VoiceHandle &voiceHandle, float pitch, float velocity, float pressure, uint16_t owner, void stealCallback(uint16_t)) {
  switch (midiType) {
    case tuningtable:
      {
        struct MidiTuningTableEntry *tte = beginTuningTableNote(pitch, velocity, pressure, owner, stealCallback);
        if (tte == nullptr) {
          return false;
        }
        voiceHandle.tte = tte;
        break;
      }
    default:
      {
        struct MpeChannelState *state = beginMpeNote(pitch, velocity, pressure, owner, stealCallback);
        if (state == nullptr) {
          return false;
        }
        state->owner = owner;
        voiceHandle.mpeState = state;
        break;
      }
  }

  voiceHandle.midiType = midiType;
  return true;
}

bool continueNote(struct VoiceHandle& voiceHandle, float pressure, uint32_t deltaUsecs) {
  switch (voiceHandle.midiType) {
    case tuningtable:
      return continueTuningTableNote(voiceHandle.tte, pressure, deltaUsecs);
    default:
      return continueMpeNote(voiceHandle.mpeState, pressure, deltaUsecs);
  }
}

bool endNote(struct VoiceHandle& voiceHandle) {
  switch (voiceHandle.midiType) {
    case tuningtable:
      return endTuningTableNote(voiceHandle.tte);
    default:
      if(endMpeNote(voiceHandle.mpeState)) {
        voiceHandle.mpeState->owner = noOne;
        return true;
      }
      return false;
  }
}

/* MIDI/MPE SETTINGS */

#define useUsbFlag          (1 << 0)  /* output MIDI on usb port */
#define useDinFlag          (1 << 1)  /* output MIDI on DIN5 port */
#define noVelocityFlag      (1 << 2)  /* don't bother sending dynamic velocity */
#define noPressureFlag      (1 << 3)  /* don't bother sending pressure */
#define skipChannel10Flag   (1 << 4)  /* some synths like to put percussion on channel 10 by default */
#define gmFlag              (1 << 5)  /* does this synth have a "general midi" mode? */
#define gm2Flag             (1 << 6)  /* does this synth have a "general midi 2" mode ? */
#define forcePbRangeFlag    (1 << 7)  /* send MIDI commands to set pitch bend range explicitly */
#define polyAfterTouchFlag  (1 << 8)  /* send poly aftertouch commands */
#define ccResetFlag         (1 << 9)  /* reset a bunch of CCs to reasonable defaults */
#define multicastTimbreFlag (1 << 10)
#define delayNoteOffFlag    (1 << 11) /* defer note-off according to pressure slew limiter */
#define noteOnFirstFlag     (1 << 12) /* send note on before pressure */
#define maxPressure126Flag  (1 << 13) /* limit max pressure to 126 */
#define bendUpOnlyFlag      (1 << 14) /* don't bend down for pitch correction, only up */
#define useETuningTableFlag (1 << 15) /* use Grey Matter E! tuning table format */ 

struct MpeSettings {
  enum midiType midiType;
  uint8_t channels;
  uint8_t pbRange; // in semitones
  uint16_t pressureBackoff;
  uint8_t pressureCC;
  uint8_t bankMsbMin;
  uint8_t bankMsbMax;
  uint8_t defaultMsb;
  uint8_t bankLsbMin;
  uint8_t bankLsbMax;
  uint8_t defaultLsb;
  uint32_t flags;
};

struct MpeSettings mpeSettingsUsbMidi = {multitimbral, 16,  2,  1000, 7,   0,  127, 0,  0,  127, 0,  useUsbFlag | polyAfterTouchFlag};
struct MpeSettings mpeSettingsXV2020  = {multitimbral, 16,  2, 15000, 7,   87, 87,  87, 64, 67,  64, useDinFlag | skipChannel10Flag | gmFlag | gm2Flag | noPressureFlag};
struct MpeSettings mpeSettingsRD300NX = {multitimbral, 16,  2,  1000, 7,   84, 121, 84, 0,  127, 0,  useDinFlag | gmFlag};
struct MpeSettings mpeSettingsFB01    = {multitimbral, 8,   2,  1000, 7,   0,  0,   0,  0,  0,   0,  useDinFlag};
struct MpeSettings mpeSettingsKSP     = {multitimbral, 4,  12,  1000, 1,   0,  0,   0,  0,  0,   0,  useDinFlag | bendUpOnlyFlag};
struct MpeSettings mpeSettingsTrinity = {multitimbral, 16,  2,  1000, 7,   0,  0,   0,  0,  3,   0,  useDinFlag | forcePbRangeFlag | gmFlag};
struct MpeSettings mpeSettingsSP300   = {multitimbral, 16,  2,  1000, 7,   0,  0,   0,  0,  0,   0,  useDinFlag};
struct MpeSettings mpeSettingsSurgeXT = {mpe,          16, 48,  1000, 128, 0,  0,   0,  0,  0,   0,  useUsbFlag | multicastTimbreFlag | noteOnFirstFlag | maxPressure126Flag };
struct MpeSettings mpeSettingsProteus = {multitimbral, 16,  2,  1000, 7,   0,  4,   4,  0,  7,   0,  useDinFlag | gmFlag | forcePbRangeFlag};
struct MpeSettings mpeSettingsMox8    = {multitimbral, 16,  2,  1400, 7,   63, 63,  63, 0,  7,   0,  useDinFlag | gmFlag | forcePbRangeFlag | ccResetFlag};
struct MpeSettings mpeSettingsPhatty  = {monovoice,    1,   2,  1000, 19,  0,  0,   0,  0,  0,   0,  useDinFlag | forcePbRangeFlag};
struct MpeSettings mpeSettingsDx7E    = {tuningtable,  1,   1,  1000, 128, 0,  0,   0,  0,  0,   0,  useDinFlag | useETuningTableFlag | noPressureFlag};
struct MpeSettings mpeSettingsDexed   = {tuningtable,  1,   1,  1000, 128, 0,  0,   0,  0,  0,   0,  useUsbFlag | noPressureFlag }

void sendMpeZones() {
  midiReadyWait();
  midiRPN(0, 6, mpeChannels, 1);
}

struct MpeSettings *currentMpeSettings = nullptr;
struct MpeSettings *mpeSettings = nullptr;

void applyMpeSettings(struct MpeSettings *settings) {
  resetAllControllers();

  useDinMidi = (settings->flags & useDinFlag) > 0;
  useUsbMidi = (settings->flags & useUsbFlag) > 0;
  midiType = settings->midiType;
  pressureBackoff = settings->pressureBackoff;
  pbRange = settings->pbRange;
  if (pressureCC != 128 && settings->pressureCC != 128) {
    ccs[pressureCC].noMulticast = false;
  }
  pressureCC = settings->pressureCC;
  if (pressureCC != 128) {
    ccs[pressureCC].noMulticast = true;
  }
  if (midiType == mpe && (settings->flags & multicastTimbreFlag ) > 0) {
    ccs[74].forceMulticast = true;
  } else {
    ccs[74].forceMulticast = false;
  }
  mpeBankMsbMin = settings->bankMsbMin;
  mpeBankMsbMax = settings->bankMsbMax;
  mpeBankLsbMin = settings->bankLsbMin;
  mpeBankLsbMax = settings->bankLsbMax;
  mpeBankMsb = settings->defaultMsb;
  mpeBankLsb = settings->defaultLsb;
  skipChannel10 = (settings->flags & skipChannel10Flag) > 0;
  doMpeDynamicVelocity = (settings->flags & noVelocityFlag) == 0;
  doMpeDynamicPressure = (settings->flags & noPressureFlag) == 0;
  doMpePolyAfterTouch = (settings->flags & polyAfterTouchFlag) != 0;
  forcePbRange = (settings->flags & forcePbRangeFlag) != 0;
  delayNoteOff = (settings->flags & delayNoteOffFlag) != 0;
  maxMpePressure = (settings->flags & maxPressure126Flag) != 0 ? 126 : 127;
  noteOnFirst = (settings->flags & maxPressure126Flag) != 0;
  bendUpOnly = (settings->flags & bendUpOnlyFlag) != 0;

  switch (settings->midiType) {
    case multitimbral:
    case monovoice:
    case tuningtable:
      firstMpeChannel = 0;
      mpeChannels = settings->channels;
      doMsb = true;
      doLsb = true;
      doProgramChange = true;
      break;
    case mpe:
      firstMpeChannel = 1;
      mpeChannels = 15;
      doMsb = false;
      doLsb = false;
      doProgramChange = false;
      sendMpeZones();
      break;
    default:
      Serial.println("unimplemented midi type");
  }

  if ((settings->flags & ccResetFlag) != 0) {
    mpeCCReset();
  }

  tuningTableType = noTuningTable;

  if ((settings->flags & useETuningTableFlag) != 0) {
    tuningTableType = eTuningTable;
    sendETuningTable();
  }

  currentMpeSettings = settings;
  mpeSettings = settings;
}

/*
 * Stop all current notes and reset MPE state.
 * Used to prevent inconsistent state when we're switching configuration settings.
 */
void mpeStop() {
  Serial.println("mpeStop");
  for (int i = 0; i < 16; i++) {
    Serial.println("mpestate " + String(i));
    struct MpeChannelState *state = &mpeState[i];
    if (state->playing) {
      midiReadyWait();
      endMpeNote(state);

      midiReadyWait();
      if (midiType == mpe) {
        midiAfterTouch(1.0, i+1);
      } else {
        midiControlChange(pressureCC, 1.0, i+1);
      }

      midiReadyWait();
      midiPitchBend(0, i+1);
      
      state->lastPitchBend = 0;
      state->lastVolume = 128;
      state->lastFilter = 128;
      state->owner = noOne;
    }
  }
  resetAllControllers();
}

void mpeSetup() {
  for (int channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
    struct MpeChannelState *state = &mpeState[channel];
    state->channel = channel;
    state->playing = false;
    state->age = 0x7fffffff;
    state->lastNote = 0;
    state->lastVolume = 128; /* deliberately out-of-range values here */
    state->lastPolyAT = 0;
    state->volume = 1.0f;
    state->lastFilter = 128;
    state->lastBendInterval = 1.0;
    state->lastPitchBend = 0;
    state->pitchBendAge = 0xffffffff;
    state->owner = noOne;
    state->lastProgramChangeSent = 128;
    state->volumeAge = 0xffffffff;
    state->lastBankLsb = 128;
    state->lastBankMsb = 128;
    state->stealCallback = nullptr;
  }

  applyMpeSettings(&mpeSettingsSurgeXT);
}

/* Used by MOX8 device preset */
void mpeCCReset(){
  for (int channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
    midiReadyWait();
    midiControlChange(   1, 127, channel+1); // mod wheel
    midiControlChange(   5,   0, channel+1); // portamento time
    midiControlChange(   7,   0, channel+1); // volume
    midiControlChange( 0xa,  64, channel+1); // pan
    midiReadyWait();
    midiControlChange(0x1f,  64, channel+1); // eg sustain level
    midiControlChange(0x40,   0, channel+1); // sustain pedal
    midiControlChange(0x41,   0, channel+1); // portamento pedal
    midiControlChange(0x42,   0, channel+1); // sostenuto pedal
    midiReadyWait();
    midiControlChange(0x47,   0, channel+1); // resonance
    midiControlChange(0x48,  64, channel+1); // eg release
    midiControlChange(0x49,  64, channel+1); // eg attack
    midiControlChange(0x4a,  64, channel+1); // filter cutoff
    midiReadyWait();
    midiControlChange(0x5b,  32, channel+1); // reverb send
    midiControlChange(0x5d,   0, channel+1); // chorus send
  }
}

/* 
 * MPE channels that aren't "owned" anymore because they've been released still
 * need to be updated to maintain proper age.
 *
 * Also, if the MIDI bus is idle and we have some channels that are due to have program
 * change messages sent, we can send them here.
 */

int idleCount = 0;

void mpeUpdate(uint32_t deltaUsecs) {
  uint16_t channelMask = 0;
  for (int channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
    struct MpeChannelState *state = &mpeState[channel];
    state->age += deltaUsecs;
    if (state->owner == noOne) {
      if (state->playing) {
        continueMpeNote(state, 0.0f, deltaUsecs);
      }
    }

    /*
    float volume = state->volume - releaseRate * ((float)deltaUsecs/1000000.0f);
    if (volume < 0.0f) {
      volume = 0.0f;
    }
    state->volume = volume; */

    if (midiIdle()) {
      idleCount += deltaUsecs;
    } else {
      idleCount = 0;
    }

    /* We don't want to flood the synth with unnecessary midi messages when we're rapidly scrolling through patches. */
    if (midiType != mpe && idleCount > 20000 * (channel+10)) {
      prepareChannel(state);
    }

    /* We've changed the CC we're using for pressure, so reset both the old one and the new one to full on. */
    if (pressureCC != pressureCCSet) {
      midiReadyWait();
      midiControlChange(pressureCCSet, 127, channel+1);
      midiControlChange(pressureCC, 127, channel+1);
      state->lastVolume = 127;
    }

    if (mpeChannels != mpeChannelsSet) {
      if (midiType == mpe && firstMpeChannel + mpeChannels > 16) {
        mpeChannels = 15;
      }
      mpeStop();
      sendMpeZones();
      mpeChannelsSet = mpeChannels;
    }
    channelMask |= 1 << channel;
  }

  doMpeMasterPitchbend(pbUp, pbDown, deltaUsecs);

  /* process any changes to current active device preset */
  if (mpeSettings != nullptr && mpeSettings != currentMpeSettings) {
    mpeStop();
    applyMpeSettings(mpeSettings);
    currentMpeSettings = mpeSettings;
  }

  /* process bank MSB/LSB unlock */
  if (unlockBankRange && !unlockBankRangeSet) {
    mpeBankMsbMax = 127;
    mpeBankMsbMin = 0;
    mpeBankLsbMax = 127;
    mpeBankLsbMin = 0;
    unlockBankRangeSet = unlockBankRange;
  }

  /* process bank MSB/LSB lock */
  if (!unlockBankRange && unlockBankRangeSet) {
    mpeBankMsbMax = currentMpeSettings->bankMsbMax;
    mpeBankMsbMin = currentMpeSettings->bankMsbMin;
    mpeBankLsbMax = currentMpeSettings->bankLsbMax;
    mpeBankLsbMin = currentMpeSettings->bankLsbMin;
    unlockBankRangeSet = unlockBankRange;
  }

  /* handle changes to pressureCC */
  if (pressureCC != pressureCCSet) {
    if (pressureCC != 128) {
      ccs[pressureCC].noMulticast = false;
    }
    pressureCCSet = pressureCC;
    if (pressureCC != 128) {
      ccs[pressureCC].noMulticast = true;
    }
  }

  /* send out any backlog of CC messages as output buffer space allows */
  bool congested = false;

  for (int cc = 0; cc < 127 && !congested; cc++) {
    if ((ccs[cc].dirty & channelMask) > 0) {
      Serial.println("dirty ccs " + String(cc) + " mask " + String(ccs[cc].dirty));
      if (midiType == mpe && !ccs[cc].forceMulticast) {
        Serial.println("mpe multicast");
        if (midiReadyLowPriority()) {
          midiControlChange(cc, ccs[cc].value, firstMpeChannel);
          ccs[cc].dirty = 0;
          Serial.println("sent multicast cc " + String(cc) + " value " + String(ccs[cc].value));
        } else {
          congested = true;
        }
      } else {
        Serial.println("non-mpe manual multicast, nomulticast " + String(ccs[cc].noMulticast));
        for (int channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels && !congested; channel++) {
          if (((1 << channel) & ccs[cc].dirty) != 0 && !ccs[cc].noMulticast) {
            if (midiReadyLowPriority()) {
              midiControlChange(cc, ccs[cc].value, channel+1);
              ccs[cc].dirty &= ~(1 << channel);
              Serial.println("sent cc " + String(cc) + " value " + String(ccs[cc].value) + " on channel " + String(channel+1));
            } else {
              congested = true;
            }
          }
        }
      }
    }    
  }
}

/* Setup, Main Loop */

int values[maxShiftRegisterBits][adcChannels] = {};
float resistances[maxShiftRegisterBits][adcChannels] = {};
float preCalibrationResistances[maxShiftRegisterBits][adcChannels] = {};
float forces[maxShiftRegisterBits][adcChannels] = {};

/*
 * Some current flows through velostat to adjacent channels as if each channel
 * were connected to the others through a resistor, whose value changes dynamically.
 * Here we store the resistance (through all paths, not just the one direct path)
 * between each channel pair.
 * The actual value stored is the reciprocal of the resistance, so we can avoid
 * floating point division later.
 */
float allPathsCalibrationMatrix[adcChannels][adcChannels] = {};
float calibrationMatrix[adcChannels][adcChannels] = {};

float lerp(float a, float b, float c){
  float range = c-a;
  float vector = b-a;
  float value = vector / range;
  if (value < 0.0) {
    return 0.0f;
  } else if (value > 1.0f) {
    return 1.0f;
  } else {
    return value;
  }
}

float lerpNoClamp(float a, float b, float c) {
  float range = c-a;
  float vector = b-a;
  return vector / range;
}

int clamp8(int a) {
  if (a < 0) {
    return 0;
  } else if (a > 255) {
    return 255;
  } else {
    return a;
  }
}

float clamp(float a) {
  if (a > 1.0f) {
    return 1.0f;
  } else if (a < 0.0f) {
    return 0.0f;
  }

  return a;
}

double clampDoubleUnit(double a) {
  if (a > 1.0) {
    return 1.0;
  } else if (a < 0.0) {
    return 0.0;
  }

  return a;
}

double clampDoubleBipolar(double a) {
  if (a > 1.0) {
    return 1.0;
  } else if (a < -1.0) {
    return -1.0;
  }

  return a;
}


enum PitchType {
  ratio,
  frequency,
  cents,
  edo
};

struct RatioPitch {
  int a;
  int b;
};

struct EdoPitch {
  int edoBase;
  int steps;
};

struct Pitch {
  enum PitchType type;
  union {
    RatioPitch ratio;
    float frequency;
    float cents;
    EdoPitch edo;
  };
};

enum KeyState {
  idle,
  playing,
  releasing,
  stolen
};

struct Key {
  int index;
  struct Pitch pitch;
  enum KeyState state;
  double intensity;
  struct VoiceHandle voiceHandle;
  float lastPressure;
};

enum knobState {
  disabled,        /* no function is assigned to this knob */
  uninitialized,   /* we have not yet read an initial value */
  inactive,        /* we're suppressing output until the user wiggles the knob the first time */
  active           /* knob is in use */
};

struct Knob {
  Knob() {};
  float max = 9000.0;
  float min = 300.0;

  float valueUpperBound = 100.0;
  float valueLowerBound = 0.0;
  float hysteresis = 500.0;
  float current;
  uint32_t data = 0;
  float initial = 0.0f;
  enum knobState state = disabled;
  bool activeByDefault = false;
  void (*action)(uint32_t data, float current) = nullptr;
  String label = "";
};

void knobCCAction(uint32_t cc, float value) {
  mpeMulticastCC(cc, value);
}

void knobReleaseRateAction(uint32_t unused, float value) {
  if (value < 0.01f) {
    value = 0.01f;
  }

  releaseRate = 0.3f / value;
}

void knobPressureExponentAction(uint32_t unused, float value) {
  value = clamp(value);

  float min = 0.25f;
  float max = 2.0f;

  pressureExponent = min + (value * (max - min));
}

#define maxKeys 113
struct Key keys[maxKeys];
int keyAllocIdx = 0;

#define maxKnobs 10
struct Knob knobs[maxKnobs];

#define controlNameLen (13)

enum ControlType {
  pressure,
  pot,
  analogIn,
  onOffSwitch,
  resistor,
  unused
};


struct Control {
  Control() {
    type = unused;
    updateFrequency = 1;
    strncpy(name, "unused", controlNameLen);
    name[controlNameLen-1] = '\0';
    bit = 0;
    channel = 0;
    updateFrequency = 0;
    key = nullptr;
    update = nullptr;
    delay = 0;
    thresholdPressure = 0;
    maxPressure = 4095;
  };

  Control(enum ControlType type, int bit, int channel, const char name_[controlNameLen], uint16_t thresholdPressure, uint16_t maxPressure):
      type{type}, bit{bit}, channel{channel}, thresholdPressure{thresholdPressure}, maxPressure{maxPressure} {
    updateFrequency = 1;
    strncpy(name, name_, controlNameLen);
    name[controlNameLen-1] = '\0';
    updateFrequency = 1;
    key = nullptr;
    update = nullptr;
    delay = 0;
    data = 0;
    held = false;
  };

  enum ControlType type;
  int bit;
  int channel;
  int updateFrequency;
  uint32_t delay;
  struct Key *key;
  char name[controlNameLen];
  void (*update)(struct Control* control, uint32_t deltaUsecs);
  uint16_t thresholdPressure;
  uint16_t maxPressure;
  uint32_t data;
  bool held;
  float area = 1.0f;
};

enum SensorType {
  sensitronics,
  velostat,
  bare
};

enum SensorType sensorType = velostat;

void stealCallback(uint16_t owner) {
  keys[owner].state = stolen;
}

void keyUpdate(struct Control* control, uint32_t deltaUsecs) {
  struct Key *key = control->key;
  float pressure = forces[control->bit][control->channel];

  float lastPressure = key->lastPressure;
  key->lastPressure = pressure;

  float velocity = mpeStaticVelocity;

  if (doMpeDynamicVelocity) {
    float delta = pressure - lastPressure;
    velocity = delta > 0
      ? pow( ( delta * 5000.0f) / (float)deltaUsecs, 0.3f)
      : pow( (-delta * 5000.0f) / (float)deltaUsecs, 0.3f);
    if (velocity > 1.0f) {
      velocity = 1.0f;
    } else if (velocity < 0.0f) {
      velocity = 0.0f;
    }
  }

  float intensity = 1.0f;
  if (doMpeDynamicPressure || doMpePolyAfterTouch) {
    if (pressure > 0.01f) {
      if (pressure < 1.0f) {
        intensity = pow((pressure * 1.01f) - 0.01f, pressureExponent);
      }
    } else {
      intensity = 0.0f;
    }
  }

  if (pressure > 0.01f && (lastPressure > 0.01f || pressure >= 1.0f || !doMpeDynamicVelocity)) {
    switch (key->state) {
      case idle:
      case releasing:
        if (!beginNote(key->voiceHandle, (double)key->pitch.ratio.a / (double)key->pitch.ratio.b, velocity, intensity, key->index, stealCallback)) {
          Serial.println("beginNote failed");
          return;
        }
        key->intensity = intensity;
        key->state = playing;
        Serial.println("noteOn " + String(key->pitch.ratio.a) + "/" + String(key->pitch.ratio.b) + " v=" + String(velocity));
        break;
      case playing:
        if (key->voiceHandle.stolen(key->index)) {
          key->state = idle;
          //key->mpeState = nullptr;
          break;
        }
        key->intensity = intensity;
        if (!continueNote(key->voiceHandle, intensity, deltaUsecs)) {
          //key->mpeState->owner = noOne;
          //key->mpeState = nullptr;
          key->state = releasing;
        }
        break;
      case stolen:
        break;
    }
  } else if (pressure <= 0.0f && (lastPressure <= 0.0f || !doMpeDynamicVelocity)) {
    switch (key->state) {
      case idle:
      case releasing:
        break;
      case playing:
        if (delayNoteOff && doMpeDynamicPressure) {
          key->intensity = intensity;
          if (!continueNote(key->voiceHandle, intensity, deltaUsecs)) {
            //key->mpeState->owner = noOne;
            //key->mpeState = nullptr;
            key->state = releasing;
          }
        } else {
          if (endNote(key->voiceHandle)) {
            //key->mpeState->owner = noOne;
            //key->mpeState = nullptr;
            key->state = releasing;
          } else {
            Serial.print("!");
          }
          key->intensity = 0.0;
        }
        break;
      case stolen:
        //key->mpeState = nullptr;
        key->state = idle;
        break;
    }
  }
}

void pbUpUpdate(struct Control* control, uint32_t deltaUsecs) {
  float force = forces[control->bit][control->channel];
  if (force < 0.0f) {
    force = 0.0f;
  }
  
  pbUp = force;
}

void pbDownUpdate(struct Control* control, uint32_t deltaUsecs) {
  float force = forces[control->bit][control->channel];
  if (force < 0.0f) {
    force = 0.0f;
  }

  pbDown = force;
}

bool debounce(struct Control* control, uint32_t deltaUsecs) {
  if (deltaUsecs >= control->delay) {
    control->delay = 0;
    return false;
  }
  
  control->delay -= deltaUsecs;
  return true;
}

#define check_debounce {if (debounce(control, deltaUsecs)) {return;}}

void menuButtonUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;
  int button = control->data;
  float force = forces[control->bit][control->channel];
  if (force > 0.02) {
    if (!control->held) {
      menuPress(button, force, deltaUsecs);
      control->held = true;
      Serial.println("menu button " + String(button + 1) + " pressed");
    }
  }
  
  if (force <= 0.0) {
    if (control->held) {
      menuRelease(button, force, deltaUsecs);
      control->held = false;
      Serial.println("menu button " + String(button + 1) + " released");
    }
  }
}

void knobUpdate(struct Control* control, uint32_t deltaUsecs) {
  int knobIndex = control->data;
  struct Knob *knob = &knobs[knobIndex];
  float r = resistances[control->bit][control->channel];
  float h = knob->hysteresis;
  bool newval = false;

  switch (knob->state) {
    case disabled:
      return;
    case uninitialized:
      if (knob->activeByDefault) {
        knob->state = active;
      } else {
        knob->initial = r;
        knob->state = inactive;
        return;
      }
      break;
    case inactive:
      {
        float excursion = abs(r - knob->initial);
        if (excursion > knob->hysteresis * 4.0f) {
          knob->state = active;
          break;
        }
      }
      return;
    case active:
      break;
  }

  if (r > knob->max) {
    r = knob->max;
  }

  if (r < knob->min) {
    r = knob->min;
  }

  if (r > knob->valueUpperBound) {
    knob->valueUpperBound = r + (h * 0.10f);
    knob->valueLowerBound = r - (h * 0.90f);
    newval = true;
  }

  if (r < knob->valueLowerBound) {
    knob->valueLowerBound = r - (h * 0.10f);
    knob->valueUpperBound = r + (h * 0.90f);
    newval = true;
  }

  if (newval) {
    float current = (r - knob->min) / (knob->max - knob->min);
    if (current > 1.0f) {
      current = 1.0f;
    } else if (current < 0.0f) {
      current = 0.0f;
    }
    knob->current = current;

    Serial.println("knob " + String(knobIndex + 1) + " value " + String(knob->current));

    if (knob->action != nullptr) {
      knob->action(knob->data, current);
    }

    status1TextUpdate(knob->label, 500000);
  }
}

void scrollUpUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  if (menuStackPos == 0) {
    return;
  }

  struct MenuItem* menu = menuStack[menuStackPos - 1];
  if (menu->childrenExtended == nullptr) {
    return;
  }

  float force = forces[control->bit][control->channel];
  if (force > 0.0f) {
    menuScroll(-1);
    control->delay = (uint32_t)((1.2f - clamp(force)) * 200000.0f);
  }
}

void scrollDownUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  if (menuStackPos == 0) {
    return;
  }

  struct MenuItem* menu = menuStack[menuStackPos - 1];
  if (menu->childrenExtended == nullptr) {
    return;
  }

  float force = forces[control->bit][control->channel];
  if (force > 0.0f) {
    menuScroll(1);
    control->delay = (uint32_t)((1.2f - clamp(force)) * 200000.0f);
  }
}

void allNotesOffFast() {
  for (int channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
    while (!midiReady()) {
      delayMicroseconds(100);
    }
    midiControlChange(123, 0, channel+1);
  }
}

void allNotesOffSlow() {
  for (int channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
    for (int note = 0; note < 128; note++) {
      while (!midiReady()) {
         delayMicroseconds(100);
      }
      midiNoteOffNoCount(note, 63, channel+1);
    }
  }

  Serial.println("all notes off");
}

void resetAllControllers() {
  for (int channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
    midiReadyWait();
    midiControlChange(121, 0, channel+1);
  }
}

void statusTextUpdate() {
  double cents = pitchToCents(transpose);
  int octave = 0;

  while (cents < -50.0) {
    octave--;
    cents += 1200.0;
  }

  while (cents > 1150.0) {
    octave ++;
    cents -= 1200.0;
  }

  int closestSemitone = (cents + 50.0)/100.0;
  String name;
  switch (closestSemitone) {
     case (0):  name = "C";  break;
     case (1):  name = "C#"; break;
     case (2):  name = "D";  break;
     case (3):  name = "D#"; break;
     case (4):  name = "E";  break;
     case (5):  name = "F";  break;
     case (6):  name = "F#"; break;
     case (7):  name = "G";  break;
     case (8):  name = "G#"; break;
     case (9):  name = "A";  break;
     case (10): name = "A#"; break;
     case (11): name = "B";  break;
     case (12): name = "C";  break;
     default:   name = "?";  break;
  }

  String type;
  switch (midiType) {
    case (monotimbral):  type = "poly"; break;
    case (tuningtable):  type = "tt";  break;
    case (multitimbral): type = "midi"; break;
    case (monovoice):    type = "mono"; break;
    case (mpe):          type = "mpe";  break;
  }

  String output = "";
  if (useUsbMidi) {
    output = output + " usb";
  }

  if (useDinMidi) {
    output = output + " din5";
  }

  windows[statusBar2].text = " " + name + String(octave+4) + " " + type + output;
  windows[statusBar2].redraw = true;
}

void status1TextUpdate(String text, uint32_t usecs) {
  windows[statusBar1].text = text;
  windows[statusBar1].textUsecs = usecs;
  windows[statusBar1].redraw = true;
}

void transposeUpdate(struct Control* control, uint32_t deltaUsecs, float transposeAmount) {
  check_debounce;

  float force = forces[control->bit][control->channel];
  if (force > 0.02f) {
    if (!control->held) {
      transpose *= transposeAmount;
      if (transpose > 8.0f) {
        transpose = 8.0f;
      } else if (transpose < 0.25f) {
        transpose = 0.25f;
      }

      statusTextUpdate();
      setPitchReference(pitchReferenceHz());
      control->held = true;

      /* re-generate and send the whole tuning table */
      if (midiType == tuningtable) {
        if (tuningTableType == eTuningTable) {
          sendETuningTable(0);
        }
      }
    }
  }
  
  if (force <= 0.0) {
    if (control->held) {
      control->held = false;
      control->delay = 200000;
    }
  }
}

void transposeUpUpdate(struct Control* control, uint32_t deltaUsecs) {
  transposeUpdate(control, deltaUsecs, 2.0f);
}

void transposeDownUpdate(struct Control* control, uint32_t deltaUsecs) {
  transposeUpdate(control, deltaUsecs, 0.5f);
}

void transposeUpSemitoneUpdate(struct Control* control, uint32_t deltaUsecs) {
  transposeUpdate(control, deltaUsecs, semitone);
}

void transposeDownSemitoneUpdate(struct Control* control, uint32_t deltaUsecs) {
  transposeUpdate(control, deltaUsecs, 1.0f/semitone);
}

void editValueIncrementUpdate(struct Control* control, uint32_t deltaUsecs) {
  if (editItem == nullptr) {
    return;
  }

  check_debounce;

  float force = forces[control->bit][control->channel];
  if (force > 0.0f) {
    auto val = *(uint32_t*)(editItem->data);
    if (editItem->maxValue == nullptr || val < *editItem->maxValue) {
      val++;
      windows[visualizerWindow].text = String(val);
      *(uint32_t*)(editItem->data) = val;
      windows[visualizerWindow].redraw = true;
      control->delay = (uint32_t)((1.0 - force) * 200000.0);
    }
  }
}

void editValueDecrementUpdate(struct Control* control, uint32_t deltaUsecs) {
  if (editItem == nullptr) {
    return;
  }

  check_debounce;

  float force = forces[control->bit][control->channel];
  if (force > 0.0f) {
    auto val = *(uint32_t*)(editItem->data);
    if (editItem->minValue == nullptr || val > *editItem->minValue) {
      val--;
      windows[visualizerWindow].text = String(val);
      *(uint32_t*)(editItem->data) = val;
      windows[visualizerWindow].redraw = true;
      control->delay = (uint32_t)((1.0 - force) * 200000.0);
    }
  }
}

struct Control controls[maxShiftRegisterBits][adcChannels];

#define CONTROL(type, bit, channel, name) (controls[bit][channel] = Control(type, bit, channel, name, thresholdPressure, maxPressure))

void controlSetupController(uint16_t thresholdPressure, uint16_t maxPressure) {
  for (int i=0; i<4; i++) {
    for (int j=0; j<adcChannels; j++) {
      controls[i][j] = Control();
    }
  }

  char const *cal;

  if (hwversion < 3) {
    cal = "cal-10k";
  } else {
    cal = "cal-2k";
  }

  CONTROL(resistor, 0, 0, cal);
  CONTROL(pressure, 1, 0, "nav-back");
  CONTROL(pressure, 2, 0, "nav-forward");
  CONTROL(pressure, 3, 0, "nav-ok");
  CONTROL(pressure, 4, 0, "nav-cancel");
  CONTROL(pressure, 5, 0, "nav-up");
  CONTROL(pressure, 6, 0, "nav-right");
  CONTROL(pressure, 7, 0, "nav-left");

  CONTROL(resistor, 0, 1, cal);
  CONTROL(pressure, 1, 1, "nav-down");
  CONTROL(pressure, 2, 1, "button-1");
  CONTROL(pressure, 3, 1, "button-2");
#if (hwversion < 2)
  CONTROL(pressure, 4, 1, "button-4");
  CONTROL(pressure, 5, 1, "button-3");
#else
  CONTROL(pressure, 4, 1, "button-3");
  CONTROL(pressure, 5, 1, "button-4");
#endif
  CONTROL(pressure, 6, 1, "button-5");
  CONTROL(pressure, 7, 1, "cal-10k");

  CONTROL(resistor, 0, 2, cal);
  CONTROL(analogIn, 1, 2, "analog-in-4");
  CONTROL(analogIn, 2, 2, "analog-in-3");
  CONTROL(analogIn, 3, 2, "analog-in-2");
  CONTROL(analogIn, 4, 2, "analog-in-1");
  CONTROL(pot, 5, 2, "rv1");
  CONTROL(pot, 6, 2, "rv2");
  CONTROL(pot, 7, 2, "rv3");

  CONTROL(resistor, 0, 3, cal);
  CONTROL(pot, 1, 3, "rv4");
  CONTROL(pot, 2, 3, "rv5");
  CONTROL(pot, 3, 3, "rv6");
  CONTROL(pot, 4, 3, "rv7");
  CONTROL(pot, 5, 3, "rv8");
  CONTROL(pot, 6, 3, "rv9");
  CONTROL(pot, 7, 3, "rv10");

  /* control button presses should be a little less sensitive */
  thresholdPressure += thresholdPressure/2;

  controls[2][1].update = menuButtonUpdate;
  controls[3][1].update = menuButtonUpdate;
  controls[4][1].update = menuButtonUpdate;
  controls[5][1].update = menuButtonUpdate;
  controls[6][1].update = menuButtonUpdate;

  controls[2][1].data = menuText1;
  controls[3][1].data = menuText2;
  controls[4][1].data = menuText3;
  controls[5][1].data = menuText4;
  controls[6][1].data = menuText5;

  controls[5][0].update = scrollUpUpdate;  // up
  controls[1][1].update = scrollDownUpdate;  // down

  controls[1][0].update = menuButtonUpdate;
  controls[1][0].data = backText;
  controls[6][0].update = editValueIncrementUpdate; // right
  controls[7][0].update = editValueDecrementUpdate; // left

  controls[5][2].update = knobUpdate;
  controls[5][2].data = 0;
  controls[6][2].update = knobUpdate;
  controls[6][2].data = 1;
  controls[7][2].update = knobUpdate;
  controls[7][2].data = 2;

  controls[1][3].update = knobUpdate;
  controls[1][3].data = 3;
  controls[2][3].update = knobUpdate;
  controls[2][3].data = 4;
  controls[3][3].update = knobUpdate;
  controls[3][3].data = 5;
  controls[4][3].update = knobUpdate;
  controls[4][3].data = 6;
  controls[5][3].update = knobUpdate;
  controls[5][3].data = 7;
  controls[6][3].update = knobUpdate;
  controls[6][3].data = 8;
  controls[7][3].update = knobUpdate;
  controls[7][3].data = 9;

  knobs[3].data = 1; /* mod wheel */
  knobs[3].action = &knobCCAction;
  knobs[3].state = active;
  knobs[3].label = "mod";

  knobs[7].data = 74; /* filter cutoff / MPE timbre */
  knobs[7].action = &knobCCAction;
  knobs[7].state = active;
  knobs[7].label = "timbre";

  knobs[8].data = 71; /* filter resonance */
  knobs[8].action = &knobCCAction;
  knobs[8].state = uninitialized;
  knobs[8].label = "resonance";

  knobs[9].data = 11; /* filter resonance */
  knobs[9].action = &knobCCAction;
  knobs[9].state = uninitialized;
  knobs[9].label = "volume";

  knobs[6].data = 91; /* reverb send */
  knobs[6].action = &knobCCAction;
  knobs[6].state = uninitialized;
  knobs[6].label = "reverb";

  knobs[0].data = 0;
  knobs[0].action = &knobReleaseRateAction;
  knobs[0].state = active;
  knobs[0].label = "slew";

  knobs[1].data = 0;
  knobs[1].action = &knobPressureExponentAction;
  knobs[1].state = active;
  knobs[1].label = "curve";
}


#define CONTROL_KEY(n, a_, b_) \
  { \
    char buf[controlNameLen]; \
    int a = a_; \
    int b = b_; \
    for(int j=0; j<oct; j++) { \
      if (b%2 == 0) { \
        b /= 2; \
      } else { \
        a *= 2; \
      } \
    } \
    snprintf(buf, controlNameLen, "%d/%d", a, b); \
    int bit = (n%8) + 8 + (oct*8); \
    int channel = (n/8) % 4; \
    CONTROL(pressure, bit, channel, buf); \
    struct Key *key = &keys[keyAllocIdx]; \
    key->index = keyAllocIdx; \
    key->pitch.type = ratio; \
    key->pitch.ratio.a = a; \
    key->pitch.ratio.b = b; \
    key->intensity = 0.0; \
    controls[bit][channel].key = key; \
    controls[bit][channel].update = keyUpdate; \
    keyAllocIdx++; \
  }

struct MidiTuningTableEntry keybedTuningTable[maxKeys];

void keybedTuningTableSetup() {
  int firstKey = middleC - (28 * 2);
  //int lastKey = middleC + (28 * 2);

  for (int i = 0; i < maxKeys; i++) {
    struct MidiTuningTableEntry *tte = &keybedTuningTable[i];
    struct Key *key = &keys[i];

    tte->noteNumber = firstKey + i;
    tte->channel = 0;
    tte->pitch = (float)key->pitch.ratio.a / (float)key->pitch.ratio.b;
    tte->users = 0;
    tte->lastPressure = 128;
    tte->pressureAge = 0x7fffffff;
  }

  midiTuningTable = keybedTuningTable;
  midiTuningTableSize = maxKeys;
}

void controlSetupKeybed(uint16_t thresholdPressure, uint16_t maxPressure) {
  int oct = 0;
  CONTROL(resistor, 8, 0, "id0-10k");
  CONTROL_KEY(1, 1, 4);
  CONTROL_KEY(2, 4, 15);
  CONTROL_KEY(3, 15, 56);
  CONTROL_KEY(4, 5, 18);
  CONTROL_KEY(5, 9, 32);
  CONTROL_KEY(6, 2, 7);
  CONTROL_KEY(7, 7, 24);

  CONTROL(resistor, 8, 1, "id1-0k");
  CONTROL_KEY(9, 3, 10);
  CONTROL_KEY(10, 5, 16);
  CONTROL_KEY(11, 9, 28);
  CONTROL_KEY(12, 21, 64);
  CONTROL_KEY(13, 1, 3);
  CONTROL_KEY(14, 27, 80);
  CONTROL_KEY(15, 7, 20);

  CONTROL(resistor, 8, 2, "id2-0k");
  CONTROL_KEY(17, 45, 128);
  CONTROL_KEY(18, 3, 8);
  CONTROL_KEY(19, 8, 21);
  CONTROL_KEY(20, 7, 18);
  CONTROL_KEY(21, 2, 5);
  CONTROL_KEY(22, 5, 12);
  CONTROL_KEY(23, 27, 64);

  CONTROL(resistor, 8, 3, "id3-0k");
  CONTROL_KEY(25, 3, 7);
  CONTROL_KEY(26, 7, 16);
  CONTROL_KEY(27, 4, 9);
  CONTROL_KEY(28, 9, 20);
  CONTROL_KEY(29, 15, 32);
  CONTROL_KEY(30, 27, 56);
  CONTROL_KEY(31, 63, 128);

  for (int oct=1; oct < 4; oct++) {
    CONTROL_KEY(2, 1, 4);
    CONTROL_KEY(3, 4, 15);
    CONTROL_KEY(4, 15, 56);
    CONTROL_KEY(5, 5, 18);
    CONTROL_KEY(6, 9, 32);
    CONTROL_KEY(7, 2, 7);
  
    CONTROL_KEY(8, 7, 24);
    CONTROL_KEY(9, 3, 10);
    CONTROL_KEY(10, 5, 16);
    CONTROL_KEY(11, 9, 28);
    CONTROL_KEY(12, 21, 64);
    CONTROL_KEY(13, 1, 3);
    CONTROL_KEY(14, 27, 80);
    CONTROL_KEY(15, 7, 20);

    CONTROL_KEY(16, 45, 128);
    CONTROL_KEY(17, 3, 8);
    CONTROL_KEY(18, 8, 21);
    CONTROL_KEY(19, 7, 18);
    CONTROL_KEY(20, 2, 5);
    CONTROL_KEY(21, 5, 12);
    CONTROL_KEY(22, 27, 64);
    CONTROL_KEY(23, 3, 7);

    CONTROL_KEY(24, 7, 16);
    CONTROL_KEY(25, 4, 9);
    CONTROL_KEY(26, 9, 20);
    CONTROL_KEY(27, 15, 32);
    CONTROL_KEY(28, 27, 56);
    CONTROL_KEY(29, 63, 128);
  }

  {
    int oct = 3;
    CONTROL_KEY(1, 1, 2);
  }
  CONTROL(pressure, 8+8+6,  3, "control-1");
  CONTROL(pressure, 8+8+7,  3, "control-2");
  CONTROL(pressure, 8+16+6, 3, "spacebar-up");
  CONTROL(pressure, 8+16+7, 3, "spacebar-dn");
  CONTROL(pressure, 8+24+6, 3, "control-4");
  CONTROL(pressure, 8+24+7, 3, "control-3");
  

  controls[8+8+6][3].update = transposeDownUpdate;
  controls[8+8+7][3].update = transposeUpUpdate;
  controls[8+24+7][3].update = transposeDownSemitoneUpdate;
  controls[8+24+6][3].update = transposeUpSemitoneUpdate;

  /* "spacebar" has huge surface area, so requires a much higher threshold */
  int pbThresholdPressure = thresholdPressure * 3.5;
  int pbMaxPressure = (thresholdPressure + maxPressure) - (pbThresholdPressure + 500);

  controls[8+16+6][3].update = pbUpUpdate;
  controls[8+16+6][3].thresholdPressure = pbThresholdPressure;
  controls[8+16+6][3].maxPressure = pbMaxPressure;
  controls[8+16+6][3].area = 4.0f;
  controls[8+16+7][3].update = pbDownUpdate;
  controls[8+16+7][3].thresholdPressure = pbThresholdPressure;
  controls[8+16+7][3].maxPressure = pbMaxPressure;
  controls[8+16+7][3].area = 4.0f;
}

void allNotesOffSlowAction(void *data) {
  allNotesOffSlow();
}

void mpeHandshakeAction(void *data) {
  if (midiType != mpe) {
    mpeStop();
    midiType = mpe;
  }

  sendMpeZones();
}

void multiTimbralAction(void *data) {
  if (midiType != multitimbral) {
    mpeStop();
    midiType = multitimbral;
  }
}

void applyMpeSettingsAction(void *data) {
  struct MpeSettings *settings = (struct MpeSettings *)data;

  mpeStop();
  applyMpeSettings(settings);
}

void sendETuningTableAction(void *data) {
  sendETuningTable(0);
}

bool enableVisualizer = true;
bool lastEnableVisualizer = enableVisualizer;
bool debugShowResistances = true;
bool debugShowCalibration = false;

struct MenuItem allNotesOffSlowMenuItem("notes off", allNotesOffSlowAction);
struct MenuItem useUsbMenuItem("USB MIDI", toggle, &useUsbMidi);
struct MenuItem useDinMenuItem("DIN5 MIDI", toggle, &useDinMidi);
struct MenuItem screen10MenuItem("10%", selection, &brightness, 25);
struct MenuItem screen25MenuItem("25%", selection, &brightness, 63);
struct MenuItem screen50MenuItem("50%", selection, &brightness, 127);
struct MenuItem screen75MenuItem("75%", selection, &brightness, 191);
struct MenuItem screen100MenuItem("100%", selection, &brightness, 255);

struct MenuItem pb2MenuItem("200 cents", selection, &pbRange, 2);
struct MenuItem pb7MenuItem("700 cents", selection, &pbRange, 7);
struct MenuItem pb12MenuItem("1200 cents", selection, &pbRange, 12);
struct MenuItem pb24MenuItem("2400 cents", selection, &pbRange, 24);
struct MenuItem pb48MenuItem("4800 cents", selection, &pbRange, 48);

struct MenuItem mpeHandshakeMenuItem("mpe init", mpeHandshakeAction);

struct MenuItem doVelocityMenuItem("velocity", toggle, &doMpeDynamicVelocity);
struct MenuItem doPressureMenuItem("pressure", toggle, &doMpeDynamicPressure);
struct MenuItem doPolyAfterTouchMenuItem("poly at", toggle, &doMpePolyAfterTouch);
struct MenuItem pressureBackoffMenuItem("p backoff", value, &pressureBackoff);

struct MenuItem mpeBankLsbMenuItem("bank LSB", value, &mpeBankLsb, &mpeBankLsbMin, &mpeBankLsbMax);
struct MenuItem mpeBankMsbMenuItem("bank MSB", value, &mpeBankMsb, &mpeBankMsbMin, &mpeBankMsbMax);

struct MenuItem debugShowResistancesMenuItem("show res", toggle, &debugShowResistances);
struct MenuItem debugShowCalibrationMenuItem("show cal", toggle, &debugShowCalibration);

struct MenuItem unlockBankRangeMenuItem("unlock", toggle, &unlockBankRange);
struct MenuItem sendETuningTableMenuItem("tx E! tt", sendETuningTableAction);

uint32_t zero = 0;
uint32_t maxMidiValue = 127;
uint32_t maxChannels = 16;

struct MenuItem pressureCCMenuItem("pressure CC", value, &pressureCC, &zero, &maxMidiValue);
struct MenuItem mpeProgramChangeMenuItem("patch", value, &programChange, &zero, &maxMidiValue);
struct MenuItem mpeChannelsMenuItem("channels", value, &mpeChannels, &zero, &maxChannels);

struct MenuItem surgeXtPresetMenuItem("Surge XT", selection, &mpeSettings, (uint32_t)&mpeSettingsSurgeXT);
struct MenuItem kspPresetMenuItem("Keystep Pro", selection, &mpeSettings, (uint32_t)&mpeSettingsKSP);
struct MenuItem rd300PresetMenuItem("RD-300NX", selection, &mpeSettings, (uint32_t)&mpeSettingsRD300NX);
struct MenuItem xv2020PresetMenuItem("XV-2020", selection, &mpeSettings, (uint32_t)&mpeSettingsXV2020);
struct MenuItem sp300PresetMenuItem("SP-300", selection, &mpeSettings, (uint32_t)&mpeSettingsSP300);
struct MenuItem trinityPresetMenuItem("Trinity", selection, &mpeSettings, (uint32_t)&mpeSettingsTrinity);
struct MenuItem fb01PresetMenuItem("FB-01", selection, &mpeSettings, (uint32_t)&mpeSettingsFB01);
struct MenuItem dx7EPresetMenuItem("DX7 with E!", selection, &mpeSettings, (uint32_t)&mpeSettingsDx7E);
struct MenuItem moxPresetMenuItem("MOX6/8", selection, &mpeSettings, (uint32_t)&mpeSettingsMox8);
struct MenuItem proteus2000PresetMenuItem("Proteus 2k", selection, &mpeSettings, (uint32_t)&mpeSettingsProteus);
struct MenuItem phattyPresetMenuItem("Slim Phatty", selection, &mpeSettings, (uint32_t)&mpeSettingsPhatty);
struct MenuItem dexedPresetMenuItem("Dexed", selection, &MpeSettings, (uint32_t)&mpeSettingsDexed);

struct MenuItem arturiaMenu("Arturia", submenu, &kspPresetMenuItem);
struct MenuItem emuMenu("E-mu", submenu, &proteus2000PresetMenuItem);
struct MenuItem korgMenu("Korg", submenu, &sp300PresetMenuItem, &trinityPresetMenuItem);
struct MenuItem moogMenu("Moog", submenu, &phattyPresetMenuItem);
struct MenuItem rolandMenu("Roland", submenu, &rd300PresetMenuItem, &xv2020PresetMenuItem);
struct MenuItem yamahaMenu("Yamaha", submenu, &dx7EPresetMenuItem, &fb01PresetMenuItem, &moxPresetMenuItem);

struct MenuItem* brandsMenu[] = {&arturiaMenu, &dexedPresetMenuItem, &emuMenu, &korgMenu, &moogMenu, &rolandMenu, &surgeXtPresetMenuItem, &yamahaMenu};
struct MenuItem outputPresetsMenu("dev presets", submenu, &brandsMenu[0], 7);
struct MenuItem pbRangeMenu("bend range", submenu, &pb2MenuItem, &pb7MenuItem, &pb12MenuItem, &pb24MenuItem, &pb48MenuItem);
struct MenuItem noteOnFirstMenuItem("note-on 1st", toggle, &noteOnFirst);
struct MenuItem maxPressureMenuItem("max p", value, &maxMpePressure, &zero, &maxMidiValue);

struct MenuItem* outputMenuItems[] = {&useUsbMenuItem, &useDinMenuItem, &outputPresetsMenu, &mpeHandshakeMenuItem, &mpeChannelsMenuItem, &pbRangeMenu, &pressureCCMenuItem, &noteOnFirstMenuItem, &maxPressureMenuItem};
struct MenuItem outputMenu("output", submenu, &outputMenuItems[0], 9);

struct MenuItem controlsMenu("controls", submenu, &doVelocityMenuItem, &doPressureMenuItem, &doPolyAfterTouchMenuItem, &pressureBackoffMenuItem);

struct MenuItem visualizerMenuItem("visualizer", toggle, &enableVisualizer);
struct MenuItem screenBrightnessMenu("brightness", submenu, &screen10MenuItem, &screen25MenuItem, &screen50MenuItem, &screen75MenuItem, &screen100MenuItem);
struct MenuItem interfaceMenu("interface", submenu, &screenBrightnessMenu, &visualizerMenuItem);
struct MenuItem patchesMenu("patches", submenu, &mpeBankMsbMenuItem,  &mpeBankLsbMenuItem, &mpeProgramChangeMenuItem, &unlockBankRangeMenuItem);

struct MenuItem debugMenu("debug", submenu, &debugShowResistancesMenuItem, &debugShowCalibrationMenuItem);

struct MenuItem configMenu("settings", submenu, &outputMenu, &controlsMenu, &interfaceMenu, &debugMenu);

struct MenuItem rootMenu("", submenu, &configMenu, &patchesMenu, &emptyMenuItem, &emptyMenuItem, &allNotesOffSlowMenuItem);

void menuSetup() {
  menuSelect(&rootMenu, 0);
}

void printWidth12(const char* str) {
  char buf[13];
  snprintf(buf, 13, "%12s", str);
  Serial.print(buf);
}

void printIntWidth4(int i) {
  char buf[20];
  snprintf(buf, 19, "%4d", i);
  Serial.print(buf);
}

void printDoubleWidth8(double d) {
  if (isnan(d)) {
    Serial.print("   -nan-");
  } else {
    char buf[20];
    snprintf(buf, 19, "%8.0f", d);
    Serial.print(buf);
  }
}

void showValues() {
  for (int bit = 0; bit < maxShiftRegisterBits; bit++) {
    Serial.print("bit ");
    printIntWidth4(bit);

    for (int channel = 0; channel < adcChannels; channel++) {
      printWidth12(controls[bit][channel].name);
      int value = values[bit][channel];
      if (value < (4095-25) || true) {
       Serial.print(" ");
       printIntWidth4(4095 - values[bit][channel]);
      } else {
        Serial.print(" ----");
      }
      Serial.print(" ");
    }
    Serial.println("");
    if ((bit + 1) % 8 == 0) {
      Serial.println("");
    }
  }
}

void showResistances() {
  for (int bit = 0; bit < maxShiftRegisterBits; bit++) {
    Serial.print("bit ");
    printIntWidth4(bit);

    for (int channel = 0; channel < adcChannels; channel++) {
      printWidth12(controls[bit][channel].name);
      Serial.print(" ");
      printIntWidth4(4095 - values[bit][channel]);
      Serial.print(" ");
      printDoubleWidth8(preCalibrationResistances[bit][channel]);
      Serial.print(" ");
      printDoubleWidth8(resistances[bit][channel]);
      Serial.print(" ");
    }
    Serial.println("");
    if ((bit + 1) % 8 == 0) {
      Serial.println("");
    }
  }
}

void setup() {
  serialSetup();
  Serial.println("begin setup");
  Serial.println("initializing LEDs...");
  ledSetup();
  Serial.println("initializing ADCs...");
  adcSetup();
  Serial.println("initializing shift registers...");
  shiftRegisterSetup();
  Serial.println("initializing screen...");
  screenSetup();
  Serial.println("initializing CAN bus...");
  canSetup();
  Serial.println("initializing MIDI...");
  midiSetup();
  Serial.println("initializing MPE...");
  mpeSetup();
  Serial.println("initializing menu...");
  menuSetup();
  statusTextUpdate();
  Serial.println("initializing audio DAC...");
  audioSetup();
  Serial.println("initializing keybed...");
  for (int i = 0; i < adcChannels; i++) {
    for (int j = 0; j < adcChannels; j++) {
      if (i==j) {
        allPathsCalibrationMatrix[i][j] = 1000000.0;
        calibrationMatrix[i][j] = 1000000.0f;
      } else {
        allPathsCalibrationMatrix[i][j] = 0.0f;
        calibrationMatrix[i][j] = 0.0f;
      }
    }
  }

  for (int i = 0; i < maxShiftRegisterBits; i++) {
    for (int j = 0; j < adcChannels; j++) {
      values[i][j] = 0.0f;
      resistances[i][j] = 1.0f;
      forces[i][j] = 0.0f;
    }
  }

  int thresholdPressure = 0;
  int maxPressure = 4096;

#if (hwversion == 0)
  thresholdPressure = 60;
  maxPressure = 500;
#elif (hwversion == 1)
  thresholdPressure = 110;
  maxPressure = 700;
#elif (hwversion == 2)
  switch (sensorType) {
    case sensitronics:
      thresholdPressure = 70;
      maxPressure = 700;
      break;
    case velostat:
      Serial.println("using velostat settings");
      thresholdPressure = 900;
      maxPressure = 3500;
      break;
    case bare:
      thresholdPressure = 100;
      maxPressure = 800;
      break;
  }
#else
  switch (sensorType) {
    case sensitronics:
      thresholdPressure = 70;
      maxPressure = 700;
      break;
    case velostat:
      Serial.println("using velostat settings");
      thresholdPressure = 650;
      maxPressure = 3000;
      break;
    case bare:
      thresholdPressure = 50;
      maxPressure = 300;
      break;
  }
#endif

  controlSetupController(thresholdPressure, maxPressure);
  controlSetupKeybed(thresholdPressure, maxPressure);
  keybedTuningTableSetup();
  delayMicroseconds(10000);
  Serial.println("end setup");
}

uint32_t usecs = 0;
uint32_t prevUsecs = 0;
uint32_t prevTimestamp = 0;
uint32_t screenRedrawAge = 1000000;
uint32_t maxRedrawAge = 10000;

void loop() {
  static int iteration = 0;
  bool verbose = (iteration % 1000 == 0);
  debugFlags = 0;

  if (verbose) {
    uint32_t timestamp = micros();
    uint32_t delta = timestamp - prevTimestamp;

    Serial.print(iteration);
    Serial.print(" ");
    Serial.print(1000.0 / ((float)delta / 1000000.0));
    Serial.print(" midi tx buffer ");
    Serial.println(Serial5.availableForWrite());
    prevTimestamp = timestamp;

    //dbgSet(adcCalibrationDebug);
  }

  uint32_t prev = 0;
  int32_t focus = 0x7ffffff; /* select a bit to focus on, and don't do any further bit shifting (for debugging) */
  int32_t curBit = -1;

  /* no bit set yet */
  calibrateADCs(verbose && debugShowCalibration, allPathsCalibrationMatrix);
  refineCalibration(allPathsCalibrationMatrix, calibrationMatrix, 4);
  if (verbose && debugShowCalibration) {
    for (int i = 0; i < adcChannels; i++) {
      Serial.print("calibration ");
      for (int j = 0; j < adcChannels; j++) {
        float r1 = allPathsCalibrationMatrix[i][j];
        float r2 = calibrationMatrix[i][j];
        Serial.print(String(1.0f/r1) + " " + String(1.0f/r2) + "  ");
      }
      Serial.println();
    }
  }

  shiftRegisterClock();
  curBit++;

  prevUsecs = usecs;
  usecs = micros();
  uint32_t delta = usecs - prevUsecs;

  while (curBit < maxShiftRegisterBits && curBit <= focus) {
    uint32_t delay = getADCDelay(values[prev], values[curBit]);
    delayMicroseconds(delay);
    int correctionIterations = curBit < 8 ? 4 : 2; /* work harder to get low error on the pot values than the keys */

    bool verboseAdc = false; // (curBit == 5 || curBit == 8) && verbose;

    readADCs(values[curBit]);
    computeResistances(verboseAdc, values[curBit], calibrationMatrix, preCalibrationResistances[curBit], resistances[curBit], correctionIterations);

    for (int channel = 0; channel < adcChannels; channel++) {
      float area = controls[curBit][channel].area;
      forces[curBit][channel] = resistanceToForce(resistances[curBit][channel], area);
      if (enableVisualizer && (iteration+curBit) % 50 == 0) {
        float r = 0.0f, g = 0.0f, b = 0.0f;

        r = 1.0f - values[curBit][channel] / 4095.0f;

        auto type = controls[curBit][channel].type;

        if (type == pressure) {
          g = forces[curBit][channel];
        } else if (type == pot) {
          b = clamp(resistances[curBit][channel] / 10500.0f);
        }

        visualizerUpdateGraph(curBit * adcChannels + channel, r, g, b);
        lastEnableVisualizer = true;
      } else {
        if (!enableVisualizer && lastEnableVisualizer) {
          windows[visualizerWindow].redraw = true;
          lastEnableVisualizer = false;
        }
      }

      auto update = controls[curBit][channel].update;
      if (update != nullptr) {
        update(&controls[curBit][channel], delta);
      }
    }
    
    prev = curBit;

    if (curBit != focus) {
      shiftRegisterClock();
      curBit++;
    } else {
      break;
    }
  }
  if (curBit != focus) {
    shiftRegisterReset(curBit);
    curBit = -1;
  } else {
    delayMicroseconds(1000);
  }

  mpeUpdate(delta);

  if (verbose  && debugShowResistances) {
    showResistances();
  }

  if (iteration % 16 == 0) { 
    int surplus = noteOnCount - noteOffCount;
    if (surplus > 10) {
      surplus = 5;
    }

    /* todo: only set if the value changes */
    leds.setPixel(5, surplus);
    leds.show();
  }

  canUpdate();

  if (dinMidi.read()){
    midiMsgsReceived++;
    int note, channel, velocity, data, cc;
    bool noteCmd = false;

    byte cmd = dinMidi.getType();
    switch(cmd) {
      case midi::NoteOn:
        noteCmd = true;
        Serial.print("got NoteOn ");
        break;
      case midi::NoteOff:
        noteCmd = true;
        Serial.print("got NoteOff ");
        break;
      case midi::PitchBend:
        Serial.println("got pitch bend");
        break;
      case midi::ControlChange:
        cc = dinMidi.getData1();
        data = dinMidi.getData2();
        Serial.println(String("got CC ") + cc + " data " + data);
        break;
      default:
        Serial.print("midi message type ");
        Serial.println(cmd);
        break;
    }

    if (noteCmd) {
      note     = dinMidi.getData1();
      velocity = dinMidi.getData2();
      channel  = dinMidi.getChannel();

      Serial.print("channel ");   Serial.print(channel);
      Serial.print(" note ");     Serial.print(note);
      Serial.print(" velocity "); Serial.println(velocity);
    }
  }

  screenRedrawAge += delta;

  menuUpdate(delta);

  if (screenRedrawAge > maxRedrawAge) {
    renderScreen(screenRedrawAge);
    screenRedrawAge = 0;
  }

  if (verbose) {
    Serial.println(String("MIDI Sent: ") + midiMsgsSent + " received: " + midiMsgsReceived + " release rate " + String(releaseRate) + " pressure exponent " + String(pressureExponent) + " pb +/- " + String(pbUp) + "  " + String(pbDown));
  }
  iteration++;
}
