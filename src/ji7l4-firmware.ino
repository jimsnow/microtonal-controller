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


#define hwversion 2

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

#define screenCSPin 255
#endif


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

/*
 * Compensate for tendency for values to "bleed over" between readings
 * on adjacent channels.
 *
 * This is just a crude heuristic, it's superseded on version 3 and later of the controller PCB
 * by calibrateADCs() and applyCalibration(). 
 */
int deaverage(int a, int b) {
  int delta = (a - b) / 10;
  if (a <= b) {
    return a;
  }
  a += delta;
  
  if (a < 0) {
    return 0;
  }
  
  if (a > 4095) {
    return 4095;
  }
  return a;
}

float adcScale = 1.0/4096.0f;

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
        cal[i][j] = 1000000.0;
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

      delayMicroseconds(100);

      float iadc = ((float)readADC(i)) / 4096.0;
      float jadc = ((float)readADC(j)) / 4096.0;

      float iv = 1.0 - (((1.0 - iadc) * 1220.0) / 1020.0);
      float jv = (jadc * 1220.0) / 1020.0;
      float rijTotal = (2440.0 / (1.0 - (iv -jv))) - 2440.0;

      if (verbose) Serial.println("calibration i:" + String(i) + " j:" + String(j) + " iv:" + String(iv) + " jv:" + String(jv) + " rijTotal:" + String(rijTotal));

      /* store the reciprocal to we don't have to do division later */
      cal[i][j] = cal[j][i] = 1.0 / rijTotal;
    }
  }

  for (int i = 0; i < adcChannels; i++) {
    int pin = adcPullupPins[i];
    pinMode(pin, OUTPUT);
    digitalWrite(pin, HIGH); /* must come after setting mode to output, otherwise ignored */
  }

  delayMicroseconds(40);
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

void applyCalibration(bool verbose, const float cal[adcChannels][adcChannels], const float vAdc[adcChannels], float r[adcChannels], const int iterations) {
  float v[adcChannels];
  float vNext[adcChannels];
  float rNext[adcChannels];
  float vOrig[adcChannels];
  float rOrigInv[adcChannels];

  for (int i = 0; i < adcChannels; i++) {
    /* approximate first guess, compensate for 200 ohm resistor */
    v[i] = 1.0 - (( max((1.0 - vAdc[i]), 0.01)  * 1220.0) / 1020.0);
    vOrig[i] = v[i];
    rOrigInv[i] = 1.0 / r[i];
  }

  if (verbose) {
    Serial.println(String(vAdc[0]) + " " + String(vAdc[1]) + " " + String(vAdc[2]) + " " + String(vAdc[3]) +
                   "|" + String(v[0]) + " " + String(v[1]) + " " + String(v[2]) + " " + String(v[3]) +
                   "|" + String(r[0]) + " " + String(r[1]) + " " + String(r[2]) + " " + String(r[3]));
  }

  float damping = 0.6; /* actual resistance is a bit higher than the calibration numbers imply, because the electrical path is shared 3 ways */

  for (int iteration = 0; iteration < iterations; iteration++) {

    //damping *= 0.8;

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

void deaverage4(int inputs[], int outputs[]) {
  for (int i = 0; i < 4; i++) {
    outputs[i] = inputs[i];
    for (int j = 0; j < 4; j++) {
      if (i==j) {
        continue;
      }

      int delta = inputs[i] - inputs[j];
      if (delta > 0) {
        outputs[i] += delta/3; 
      }
    }

    if (outputs[i] > 4095) {
      outputs[i] = 4095;
    }
  }
}

int deblur(int prev, int cur) {
  int delta = cur - prev;
  if (delta < 0) {
    return cur;
  }

  int adjust = delta / 20;
  if (adjust > 100) {
    adjust = 100;
  }

  int updated = cur + adjust;
  if (updated > 4095) {
    return 4095;
  }

  return updated;
}

int readADC(int channel) {
  return adc->adc0->analogRead(adcPins[channel]);
}

void readADCs(bool verbose, int *values, float *resistances, const float calibrationMatrix[adcChannels][adcChannels], int iterations) {
  adc->adc0->startSingleRead(adcPins[0]);
  adc->adc1->startSingleRead(adcPins[1]);

  while(!adc->adc0->isComplete()) {};
  int value0 = adc->adc0->readSingle();
  while(!adc->adc1->isComplete()) {};
  int value1 = adc->adc1->readSingle();

  delayMicroseconds(1);

  adc->adc0->startSingleRead(adcPins[2]);
  adc->adc1->startSingleRead(adcPins[3]);

  while(!adc->adc0->isComplete()) {};
  int value2 = adc->adc0->readSingle();
  while(!adc->adc1->isComplete()) {};
  int value3 = adc->adc1->readSingle();

  if (hwversion < 3) {
    int preAverage [] = {value0, value1, value2, value3};
    deaverage4(preAverage, values);

    if (resistances != nullptr) {
      resistances[0] = valueToResistance(value0);
      resistances[1] = valueToResistance(value1);
      resistances[2] = valueToResistance(value2);
      resistances[3] = valueToResistance(value3);
    }
  } else {
    if (resistances != nullptr) {
      
      /* compensate for effects of conduction between channels across velostat */
      float vAdc[adcChannels] = {value0 * adcScale, value1 * adcScale, value2 * adcScale, value3 * adcScale};
      resistances[0] = valueToResistance(value0);
      resistances[1] = valueToResistance(value1);
      resistances[2] = valueToResistance(value2);
      resistances[3] = valueToResistance(value3);
      applyCalibration(verbose, calibrationMatrix, vAdc, resistances, iterations);
    }
  }

  values[0] = value0;
  values[1] = value1;
  values[2] = value2;
  values[3] = value3;
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
};

#define numWindows 11

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
  visualizerWindow,
  statusBar
};


struct Window windows[numWindows];

struct Point cursor(const struct Window &window) {
  uint16_t left = window.extent.p1.x;
  //uint16_t right = window.extent.p2.x;
  uint16_t top = window.extent.p1.y;
  uint16_t bottom = window.extent.p2.y;
  uint16_t bheight = bottom-top;

  return Point(left+6, top+(bheight/2) - 9);
}

void setWindowCursor(const struct Window &window) {
  Point c = cursor(window);
  tft.setCursor(c.x, c.y);
}

uint32_t brightness = 127;
uint32_t brightnessSet = 255;

void renderScreen() {
  //tft.fillScreen(ILI9341_BLACK);

  for (int i=0; i<numWindows; i++) {
    if (!windows[i].redraw) {
      continue;
    }

    Rectangle r = windows[i].extent;
    if (!windows[i].enabled) {
      tft.fillRect(r.p1.x+1, r.p1.y+1, (r.p2.x-r.p1.x)-2, (r.p2.y-r.p1.y)-2, ILI9341_BLACK);
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

  if (brightness != brightnessSet) {
    analogWrite(backlightPin, brightness);
    brightnessSet = brightness;
    Serial.println("set brightness to " + String(brightness) + "/255");
  }
}

void screenSetup() {
  for (int i=0; i<numWindows; i++) {
    windows[i].extent = Rectangle(0,0,0,0);
    windows[i].text = "";
    windows[i].bgcolor = ILI9341_DARKGREY;
    windows[i].fgcolor = ILI9341_WHITE;
    windows[i].enabled = true;
    windows[i].redraw = true;
    windows[i].highlight = false;
  }


  for (int i=0; i<5; i++) {
    windows[i].extent = Rectangle (0, i*menuItemHeight, menuWidth, (i+1)*menuItemHeight);
  }

  windows[5].extent  = Rectangle (menuWidth, 0, (menuWidth+navButtonWidth), menuItemHeight);
  windows[6].extent  = Rectangle (menuWidth + navButtonWidth, 0, width, menuItemHeight);
  windows[7].extent  = Rectangle (menuWidth, menuItemHeight, (menuWidth+navButtonWidth), menuItemHeight*2);
  windows[8].extent  = Rectangle (menuWidth + navButtonWidth, menuItemHeight, width, menuItemHeight*2);
  windows[statusBar].extent  = Rectangle (menuWidth, height-statusHeight, width, height);
  windows[statusBar].enabled = true;
  windows[visualizerWindow].extent = Rectangle (menuWidth, menuItemHeight*2, width, height-statusHeight);
  windows[visualizerWindow].bgcolor = ILI9341_BLACK;

  windows[backText].text = "back";
  windows[backText].enabled = false;
  windows[fwdText].text = "forward";
  windows[fwdText].enabled = false;
  windows[cancelText].text = "cancel";
  windows[cancelText].enabled = false;
  windows[okText].text = "ok";
  windows[okText].enabled = false;

  pinMode(backlightPin, OUTPUT);
  analogWriteFrequency(backlightPin, 3611*2); /* default is 3.611 kHz */
  analogWrite(backlightPin, brightness);
  brightnessSet = brightness;

  tft.begin();
  tft.setRotation(1);
  tft.setClock(100000000);
  tft.fillScreen(ILI9341_BLACK);

  renderScreen();

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
  MenuItem(String text, void (*select)(void *data), void *data = nullptr) : text{text}, type{action}, select{select}, data{data} {}
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
  uint8_t numChildren = 0;
  struct MenuItem* children[5];
  bool highlight = false;
};

struct MenuItem emptyMenuItem = MenuItem("", empty);

struct MenuItem* menu[5] = {&emptyMenuItem, &emptyMenuItem, &emptyMenuItem, &emptyMenuItem, &emptyMenuItem};

struct MenuItem* menuStack[10];
uint16_t menuStackPos = 0;

void statusTextUpdate();

void menuSelect(struct MenuItem *item, uint16_t button) {
  if (item->type == submenu) {
    for (int i=0; i < 5; i++) {
      if (i < item->numChildren) {
        windows[i].bgcolor = ILI9341_DARKGREY;
        windows[i].text = item->children[i]->text;
        windows[i].enabled = true;
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
  }

  if (item->type == toggle){
    if (item->data != nullptr) {
      *((bool*)(item->data)) = !*(bool*)(item->data);
      Serial.println("toggled value is " + String(*(bool*)(item->data)));
    } else {
      Serial.println("no value pointer to toggle");
    }
  }

  if (item->type == selection) {
    if (item->data != nullptr) {
      *((uint32_t*)(item->data)) = item->defaultValue;
    }
  }

  if (item->select != nullptr) {
    item->select(item->data);
  }

  if (menuStackPos > 0) {
    struct MenuItem *parent = menuStack[menuStackPos-1];
    for (int i = 0; i < parent->numChildren; i++) {
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
    windows[statusBar].redraw = true;
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

MIDI_CREATE_INSTANCE(HardwareSerial, Serial5, dinMidi);

uint64_t midiMsgsSent = 0;
uint64_t midiMsgsReceived = 0;

#define doMidi(func, ...) {midiMsgsSent++; if (useUsbMidi) {usbMIDI.func(__VA_ARGS__);} if (useDinMidi) {dinMidi.func(__VA_ARGS__);}}

void midiNoteOn(uint8_t note, uint8_t velocity, uint8_t channel) {
  doMidi(sendNoteOn, note, velocity, channel);
}

void midiNoteOff(uint8_t note, uint8_t velocity, uint8_t channel) {
  doMidi(sendNoteOff, note, velocity, channel);
}

void midiPitchBend(int16_t pb, uint8_t channel) {
  doMidi(sendPitchBend, pb, channel);
}

void midiAfterTouch(uint8_t volume, uint8_t channel) {
  doMidi(sendAfterTouch, volume, channel);
}

void midiControlChange(uint8_t cc, uint8_t value, uint8_t channel) {
  doMidi(sendControlChange, cc, value, channel);
}

void midiProgramChange(uint8_t bank, uint8_t channel) {
  doMidi(sendProgramChange, bank, channel);
}

int midiBufferSize = 0;

void midiSetup(){
  Serial.println("serial fifo size " + String(Serial5.availableForWrite()));
  dinMidi.begin();
  midiBufferSize = Serial5.availableForWrite();
  Serial.println("midiBufferSize set to " + String(midiBufferSize));
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

int mpeChannels = 16;
int firstMpeChannel = 0;

#define noOne 0xffff

struct MpeChannelState{
  int channel;
  bool playing;
  uint32_t age;
  uint8_t lastNote;
  uint8_t lastVolume;
  uint8_t lastFilter;
  uint16_t lastPitchBend;
  double lastBendInterval;
  uint32_t pitchBendAge;
  uint16_t owner;
  int lastProgramChangeSent;
  uint32_t volumeAge;
  uint8_t lastBankMsb;
  uint8_t lastBankLsb;
  double originalPitch;
  void (*stealCallback)(uint16_t owner);
};

struct MpeChannelState mpeState[16];
uint8_t programChange = 0;

int noteOnCount = 0;
int noteOffCount = 0;

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
    noteOffCount++;
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

#define middleC 60 /* midi note */
double pbRange = 2.0;
uint8_t pressureCC = 0;
uint8_t pressureMax = 127;

/* seconds from full to none, or vice versa */
double attack = 0.0;
double decay = 1.0;

double pitchToCents(double pitch) {
  return (log(pitch) / log(2.0)) * 1200.0;
}

double centsToPitch(double cents) {
  const double factor = pow(2.0, (1.0/1200.0));
  return pow(factor, cents);
}

double semitone = pow(2.0, (1.0/12.0));

enum midiType {
  monotimbral,
  mts,
  multitimbral,
  mpe
};

double transpose = 1.0;
int mpeBankLsbMin = 0;
int mpeBankLsbMax = 0;
int mpeBankMsbMin = 0;
int mpeBankMsbMax = 0;
int mpeBankMsb = 0;
int mpeBankLsb = 0;
enum midiType midiType = multitimbral;
int pressureBackoff = 5000;

double masterPbUpRange = 7.0/6.0;
double masterPbDownRange = 6.0/7.0;

double pbUp = 0.0;
double pbDown = 0.0;

uint32_t masterPbAge = 0;
int16_t lastMasterPb = 0;
uint32_t masterPbBackoff = 8000;

double masterPbRange = 12.0;

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

struct MpeChannelState *beginMpeNote(double pitch, int owner, double velocity, double pressure, void stealCallback(uint16_t)) {
  double shift = pitchToCents(pitch * transpose);
  //state->originalPitch = shift;

  int semitones = shift / 100.0;
  double cents = shift - (semitones * 100.0);
  if (cents > 50.0) {
    semitones++;
    cents -= 100.0;
  } else if (cents < -50.0) {
    semitones--;
    cents += 100.0;
  }

  int note = middleC + semitones;
  int16_t pb;
  if (midiType == mpe) {
    pb = calculatePitchBend(0.0, 0.0, centsToPitch(cents), pbRange);
  } else {
    pb = calculatePitchBend(pbUp, pbDown, centsToPitch(cents), pbRange);
  }

  int v = (int) (1.0 + (velocity*126));

  if (v > 127) {
    v = 127;
  } else if (v < 1) {
    v = 1;
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
    state->lastPitchBend = pb;

    if (midiType != mpe) {

      if (state->lastBankMsb != mpeBankMsb) {
        midiControlChange(0, mpeBankMsb, state->channel+1);
        state->lastBankMsb = mpeBankMsb;
      }

      if (state->lastBankLsb != mpeBankLsb) {
        midiControlChange(32, mpeBankLsb, state->channel+1);
        state->lastBankLsb = mpeBankLsb;
      }

      if (state->lastProgramChangeSent != programChange) {
        midiProgramChange(programChange, state->channel+1);
        state->lastProgramChangeSent = programChange;
        Serial.print("sent program change on channel ");
        Serial.println(state->channel);

        // set pitch bend range
        midiControlChange(100, 0, state->channel+1);
        midiControlChange(101, 0, state->channel+1);
        midiControlChange(6, (int)pbRange, state->channel+1);
        midiControlChange(38, 0, state->channel+1);
        midiControlChange(101, 127, state->channel+1);
        midiControlChange(100, 127, state->channel+1);

        /* turn down resonance */
        // midiControlChange(71, 0, state->channel+1);
        /* turn up the reverb */
        // midiControlChange(91, 127, state->channel+1);
      }
    }
    midiPitchBend(pb, state->channel+1);
    state->volumeAge = pressureBackoff;
    continueMpeNote(state, pressure, 0);

    midiNoteOn(note, v, state->channel+1);
    noteOnCount++;
    Serial.print("sent note-on on channel ");
    Serial.print(state->channel+1);
    Serial.print(" note ");
    Serial.print(note);
    Serial.print(" cents ");
    Serial.print(cents);
    Serial.print(" (");
    Serial.print(centsToPitch(cents));
    Serial.print(") pb ");
    Serial.print(pb);
    Serial.print(" range ");
    Serial.println(pbRange);
  }

  state->age = 0;

  return state;
}

void continueMpeNote(struct MpeChannelState *state, double pressure, uint32_t deltaUsecs) {

  /* rate limit pressure updates */
  if (deltaUsecs < pressureBackoff - state->volumeAge) {
    state->volumeAge += deltaUsecs;
    return;
  }

  int max = 127;

  uint8_t volume = pressure * (double)max;
  if (volume > max) {
    volume = max;
  }

  if (volume == state->lastVolume) {
    state->volumeAge = 0;
    return;
  }

  if(midiReadyLowPriority()) {
    if (midiType == mpe) {
      midiAfterTouch(volume, state->channel+1);
    } else {
      midiControlChange(pressureCC, volume, state->channel+1);
    }
    state->lastVolume = volume;
    state->volumeAge = 0;
  } else {
    Serial.print(".");
    state->volumeAge += deltaUsecs;
  }
}

bool endMpeNote(struct MpeChannelState *state) {
  if (midiReady()) {
    state->age = 0;
    /* if note wasn't already stolen */
    if (state->playing == true) {
      midiNoteOff(state->lastNote, 5, state->channel+1);
      noteOffCount++;
    }
    state->playing = false;
    Serial.print("sent note-off channel ");
    Serial.print(state->channel+1);
    Serial.print(" note ");
    Serial.println(state->lastNote);

    return true;
  }
  return false;
}

void doMpeMasterPitchbend(double pbUp, double pbDown, uint32_t deltaUsecs) {
  if (deltaUsecs < masterPbBackoff - masterPbAge) {
    masterPbAge += deltaUsecs;
    return;
  }

  if (midiType == mpe) {
    int16_t pbi = calculatePitchBend(pbUp, pbDown, 1.0, masterPbRange);
    if (pbi != lastMasterPb) {
      midiPitchBend(pbi, 1);
      lastMasterPb = pbi;
    }
  } else if (midiType == multitimbral || midiType == monotimbral) {
    for (int channel = firstMpeChannel; channel < firstMpeChannel+mpeChannels; channel++) {
      if (mpeState[channel].playing) {
        int pbi = calculatePitchBend(pbUp, pbDown, mpeState[channel].lastBendInterval, pbRange);
        if (mpeState[channel].lastPitchBend != pbi) {
          midiPitchBend(pbi, channel+1);
          mpeState[channel].lastPitchBend = pbi;
        }
      }
    }
  }

  masterPbAge = 0;
}

struct MpeSettings {
  enum midiType midiType;
  const char* name;
  int channels;
  double pbRange; // in semitones
  int pressureBackoff;
  uint8_t pressureCC;
  uint8_t bankMsbMin;
  uint8_t bankMsbMax;
  uint8_t bankLsbMin;
  uint8_t bankLsbMax;
  bool useDinMidi;
  bool skipChannel10;
};

struct MpeSettings mpeSettingsDefault = {multitimbral, "default",    16,  2.0,  5000, 7,  0,  0,  0,  0,   true, false};
struct MpeSettings mpeSettingsXV2020  = {multitimbral, "xv-2020",    16,  2.0, 20000, 7,  87, 87, 64, 255, true,  true};
struct MpeSettings mpeSettingsFB01    = {multitimbral, "fb-01",       8,  2.0,  5000, 7,  0,  0,  0,  0,   true, false};
struct MpeSettings mpeSettingsKSP     = {multitimbral, "keystep pro", 4, 12.0,  5000, 1,  0,  0,  0,  0,   true, false};
struct MpeSettings mpeSettingsTrinity = {multitimbral, "trinity",    16,  2.0,  5000, 74, 0,  0,  0,  3,   true, false};
struct MpeSettings mpeSettingsSurgeXT = {mpe,          "surge-xt",   16, 48.0,  5000, 74, 0,  0,  0,  0,  false, false};

void sendMpeZones(){
  midiReadyWait();
  midiControlChange(64, 06, 1);
  midiReadyWait();
  midiControlChange(65, 00, 1);
  midiReadyWait();
  midiControlChange(06, 15, 1);
}

struct MpeSettings *currentMpeSettings = nullptr;
struct MpeSettings *mpeSettings = nullptr;

void applyMpeSettings(struct MpeSettings *settings) {
  useDinMidi = settings->useDinMidi;
  useUsbMidi = !settings->useDinMidi;
  midiType = settings->midiType;
  pressureBackoff = settings->pressureBackoff;
  pbRange = settings->pbRange;
  pressureCC = settings->pressureCC;
  mpeBankMsbMin = settings->bankMsbMin;
  mpeBankMsbMax = settings->bankMsbMax;
  mpeBankLsbMin = settings->bankLsbMin;
  mpeBankLsbMax = settings->bankLsbMax;
  mpeBankMsb = mpeBankMsbMin;
  mpeBankLsb = mpeBankLsbMin;
  skipChannel10 = settings->skipChannel10;

  switch (settings->midiType) {
    case multitimbral:
      firstMpeChannel = 0;
      mpeChannels = settings->channels;
      break;
    case mpe:
      firstMpeChannel = 1;
      mpeChannels = 15;

      sendMpeZones();
      break;
    default:
      Serial.println("unimplemented midi type");
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
      state->lastVolume = 127;
      state->lastFilter = 127;
      state->owner = noOne;
    }
  }
  Serial.println("mpeStop resetting controllers");
  Serial.println("mpeStop midiBufferSize " + String(midiBufferSize));
  resetAllControllers();
  Serial.println("done");
}

void mpeSetup() {
  for (int channel = firstMpeChannel; channel < mpeChannels; channel++) {
    struct MpeChannelState *state = &mpeState[channel];
    state->channel = channel;
    state->playing = false;
    state->age = 0xffffffff;
    state->lastNote = 0;
    state->lastVolume = 127;
    state->lastFilter = 127;
    state->lastBendInterval = 1.0;
    state->lastPitchBend = 0;
    state->pitchBendAge = 0xffffffff;
    state->owner = noOne;
    state->lastProgramChangeSent = -1;
    state->volumeAge = 0xffffffff;
    state->lastBankLsb = 0;
    state->lastBankMsb = 0;
    state->stealCallback = nullptr;
  }

  Serial.println("mpeSetup midiBufferSize " + String(midiBufferSize));

  applyMpeSettings(&mpeSettingsSurgeXT);
  //applyMpeSettings(&mpeSettingsKSP);
  //applyMpeSettings(&mpeSettingsXV2020);
}

/* 
 * MPE channels that aren't "owned" anymore because they've been released still
 * need to be updated to maintain proper age and to apply appropriate pressure slew.
 */
void mpeUpdate(uint32_t deltaUsecs) {
  for (int channel = firstMpeChannel; channel < firstMpeChannel+mpeChannels; channel++) {
    struct MpeChannelState *state = &mpeState[channel];
    if (state->owner == noOne) {
      state->age += deltaUsecs;
    }
  }

  if (mpeSettings != nullptr && mpeSettings != currentMpeSettings) {
    Serial.println("stopping in progress notes...");
    mpeStop();
    Serial.println("changing device output settings...");
    applyMpeSettings(mpeSettings);
    Serial.println("done");
    currentMpeSettings = mpeSettings;
  }
}

/* Setup, Main Loop */

int values[maxShiftRegisterBits][adcChannels] = {};
float resistances[maxShiftRegisterBits][adcChannels] = {};

/*
 * Some current flows through velostat to adjacent channels as if each channel
 * were connected to the others through a resistor, whose value changes dynamically.
 * Here we store the resistance (through all paths, not just the one direct path)
 * between each channel pair.
 * The actual value stored is the reciprocal of the resistance, so we can avoid
 * floating point division later.
 */
float calibrationMatrix[adcChannels][adcChannels] = {};

void setLed(int led, int bit, int channel) {
  int value = values[bit][channel];

  int pressure = ((4096 - value) / 4)-25;
  if (pressure <= 0)
    pressure = 0;
  else
    pressure--;

  if (pressure > 255) pressure = 255;

  int color = (pressure << 16) | (pressure << 8) | pressure;

  leds.setPixel(led, color);
}

float lerp(int a, int b, int c){
  float range = c-a;
  float vector = b-a;
  float value = vector / range;
  if (value < 0.0) {
    return 0;
  } else if (value > 1.0) {
    return 1.0;
  } else {
    return value;
  }
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
  struct MpeChannelState *mpeState;
  int lastPressure;
};

#define maxKeys 113
struct Key keys[maxKeys];
int keyAllocIdx = 0;

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
};

enum SensorType {
  sensitronics,
  velostat,
  bare
};

enum SensorType sensorType = velostat;

// 1/seconds to fall from full intensity
double releaseRate = 2.0;

void stealCallback(uint16_t owner) {
  keys[owner].state = stolen;
}

bool doMpeDynamicVelocity = true;
bool doMpeDynamicPressure = true;
float mpeStaticVelocity = 0.75;
float minVelocity = 0.2;

void keyUpdate(struct Control* control, uint32_t deltaUsecs) {
  struct Key *key = control->key;
  int value = values[control->bit][control->channel];
  uint16_t threshold = control->thresholdPressure;
  uint16_t maxPressure = control->maxPressure;

  /* debounce */
  if (key->state == playing) {
    threshold -= 6;
  } else {
    threshold += 6;
  }

  if (key->mpeState != nullptr) {
    key->mpeState->age += deltaUsecs;
  }
  /*
  if (deltaUsecs >= control->delay) {
    control->delay = 0;
  } else {
    control->delay -= deltaUsecs;
  }*/

  int pressure = (4095-value) - threshold;
  if (pressure < 0) {
    pressure = 0;
  }

  if (pressure > maxPressure) {
    pressure = maxPressure;
  }

  int lastPressure = key->lastPressure;
  key->lastPressure = pressure;

  double delta = (pressure-lastPressure) / (double)maxPressure;
  double velocity = mpeStaticVelocity;

  if (doMpeDynamicVelocity) {
    velocity = delta > 0
      ? pow( ( (double)delta * 20000.0) / (double)deltaUsecs, 0.8)
      : pow( (-(double)delta * 20000.0) / (double)deltaUsecs, 0.8);
    if (velocity > 1.0) {
      velocity = 1.0;
    } else if (velocity < minVelocity) {
      velocity = minVelocity;
    }
  }

  double intensity = 1.0;
  if (doMpeDynamicPressure) {
    intensity = pow(pressure / (double)maxPressure, 0.4);
    double minIntensity = key->intensity - releaseRate * ((double)deltaUsecs/1000000.0);

    if (intensity < minIntensity) {
      intensity = minIntensity;
    }

    if (intensity > 1.0) {
      intensity = 1.0;
    } else if (intensity < 0.0) {
      intensity = 0.0;
    }
  }

  if (pressure > 0 && (lastPressure > 0 || velocity == 1.0)) {
    switch (key->state) {
      case idle:
      case releasing:
        Serial.print("attempting noteOn ");
        Serial.print(key->pitch.ratio.a);
        Serial.print("/");
        Serial.println(key->pitch.ratio.b);
        Serial.print(" ");
        Serial.print(pressure);
        Serial.print("-");
        Serial.print(lastPressure);
        Serial.print(" ");
        Serial.print(" velocity ");
        Serial.print(velocity);
        key->mpeState = beginMpeNote((double)key->pitch.ratio.a / (double)key->pitch.ratio.b, key->index, velocity, intensity, stealCallback);
        if (key->mpeState == nullptr) {
          return;
        }
        key->intensity = intensity;
        key->state = playing;
        break;
      case playing:
        if (key->mpeState->owner != key->index) {
          key->state = idle;
          key->mpeState = nullptr;
          break;
        }
        continueMpeNote(key->mpeState, intensity, deltaUsecs);
        key->intensity = intensity;
        break;
      case stolen:
        break;
    }
  } else if (pressure <= 0 && lastPressure <= 0) {
    switch (key->state) {
      case idle:
      case releasing:
        break;
      case playing:
      case stolen:
        if (endMpeNote(key->mpeState)) {
          key->mpeState->owner = noOne;
          key->mpeState = nullptr;
          key->state = releasing;
        } else {
          Serial.print("!");
        }
        key->intensity = 0.0;
        break;
    }
  }
    /*
    Serial.print("note held ");
    Serial.print(control->key->pitch.ratio.a);
    Serial.print("/");
    Serial.print(control->key->pitch.ratio.b);
    Serial.print(" pressure ");
    Serial.println(pressure); */
  
}

void pbUpUpdate(struct Control* control, uint32_t deltaUsecs) {
  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - control->thresholdPressure;
  if (pressure < 0) {
    pressure = 0;
  }
  //Serial.println(pressure);
  //Serial.println(control->maxPressure);
  pbUp = pow((double)pressure/(double)control->maxPressure, 0.4);
  //doMpeMasterPitchbend(pbUp, pbDown, deltaUsecs/2);
}

void pbDownUpdate(struct Control* control, uint32_t deltaUsecs) {
  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - control->thresholdPressure;
  if (pressure < 0) {
    pressure = 0;
  }
  pbDown = pow((double)pressure/(double)control->maxPressure, 0.4);
  doMpeMasterPitchbend(pbUp, pbDown, deltaUsecs);
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
  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - control->thresholdPressure;
  if (pressure > 0) {
    if (!control->held) {
      menuPress(button, pressure, deltaUsecs);
      control->held = true;
      Serial.println("menu button " + String(button) + " pressed " + String(midiBufferSize));
    }
  } else {
    if (control->held) {
      Serial.println("menu button " + String(button) + " released " + String(midiBufferSize));
      menuRelease(button, pressure, deltaUsecs);
      control->held = false;
      Serial.println("menu Release called" + String(midiBufferSize));
    }
  }
}

void incProgramChangeUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  if (programChange == 127) {
    return;
  }

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - control->thresholdPressure;
  if (pressure > 0) {
    programChange++;
    control->delay = 100000;
    Serial.print("programChange ");
    Serial.println(programChange);
  }
}

void decProgramChangeUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  if (programChange == 0) {
    return;
  }

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - control->thresholdPressure;
  if (pressure > 0) {
    programChange--;
    control->delay = 100000;
    Serial.print("programChange ");
    Serial.println(programChange);
  }
}

void allNotesOffSlowUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - control->thresholdPressure;
  if (pressure > 0) {
    for (int channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
      for (int note = 0; note < 128; note++) {
        while (!midiReady()) {
          delayMicroseconds(100);
        }
        midiNoteOff(note, 63, channel+1);
      }
    }

    control->delay = 100000;
    Serial.println("all notes off (the slow way)");
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
      midiNoteOff(note, 63, channel+1);
    }
  }

  Serial.println("all notes off");
}

void resetAllControllers() {
  for (int channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
    Serial.println("resetting channel " + String(channel+1) + " midiBufferSize " + String(midiBufferSize));
    midiReadyWait();
    Serial.println("sending control change");
    midiControlChange(121, 0, channel+1);
  }
}

void allNotesOffUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - control->thresholdPressure;
  if (pressure > 0) {
    allNotesOffSlow();
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
     case (0): name = "C"; break;
     case (1): name = "C#"; break;
     case (2): name = "D"; break;
     case (3): name = "D#"; break;
     case (4): name = "E"; break;
     case (5): name = "F"; break;
     case (6): name = "F#"; break;
     case (7): name = "G"; break;
     case (8): name = "G#"; break;
     case (9): name = "A"; break;
     case (10): name = "A#"; break;
     case (11): name = "B"; break;
     case (12): name = "C"; break;
     default: name = "?"; break;
  }

  String type;
  switch (midiType) {
    case (monotimbral): type = "mono"; break;
    case (mts): type = "mts"; break;
    case (multitimbral): type= "midi"; break;
    case (mpe): type = "mpe"; break;
  }

  String output = "";
  if (useUsbMidi) {
    output = output + " usb";
  }

  if (useDinMidi) {
    output = output + " din5";
  }

  windows[statusBar].text = " " + name + String(octave+4) + " " + type + output;
  windows[statusBar].redraw = true;
}

void transposeUpUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - control->thresholdPressure;
  if (pressure > 50) {
    transpose *= 2.0;
    if (transpose > 8.0) {
      transpose = 8.0;
    }
    control->delay = 200000;
    statusTextUpdate();
    setPitchReference(pitchReferenceHz());
  }
}

void transposeDownUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - control->thresholdPressure;
  if (pressure > 50) {
    transpose *= 0.5;
    if (transpose < 0.25) {
      transpose = 0.25;
    }
    control->delay = 200000;
    statusTextUpdate();
    setPitchReference(pitchReferenceHz());
  }
}

void transposeUpSemitoneUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - control->thresholdPressure;
  if (pressure > 50) {
    transpose *= semitone;
    if (transpose > 8.0) {
      transpose = 8.0;
    }
    control->delay = 200000;
    statusTextUpdate();
    setPitchReference(pitchReferenceHz());
  }
}

void transposeDownSemitoneUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - control->thresholdPressure;
  if (pressure > 50) {
    transpose *= 1.0/semitone;
    if (transpose < 0.25) {
      transpose = 0.25;
    }
    control->delay = 200000;
    statusTextUpdate();
    setPitchReference(pitchReferenceHz());
  }
}


void bankLsbUpUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - control->thresholdPressure;
  if (pressure > 0) {
    mpeBankLsb++;
    if (mpeBankLsb > mpeBankLsbMax) {
      mpeBankLsb = mpeBankLsbMax;
    }
    Serial.print("bank lsb ");
    Serial.println(mpeBankLsb);
    control->delay = 100000;
  }
}

void bankLsbDownUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - control->thresholdPressure;
  if (pressure > 0) {
    mpeBankLsb--;
    if (mpeBankLsb < mpeBankLsbMin) {
      mpeBankLsb = mpeBankLsbMin;
    }
    Serial.print("bank lsb ");
    Serial.println(mpeBankLsb);
    control->delay = 100000;
  }
}

void mpeZoneUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - control->thresholdPressure;
  if (pressure > 0) {
    sendMpeZones();
    Serial.println("registering MPE zones");
    control->delay = 100000;
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

  //controls[2][1].update = incProgramChangeUpdate;
  //controls[3][1].update = decProgramChangeUpdate;
  //controls[6][1].update = allNotesOffUpdate;
  //controls[4][1].update = allNotesOffSlowUpdate;

  controls[5][0].update = incProgramChangeUpdate;
  controls[1][1].update = decProgramChangeUpdate;

  controls[1][0].update = menuButtonUpdate;
  controls[1][0].data = backText;
  controls[4][0].update = bankLsbUpUpdate;
  //controls[5][1].update = mpeZoneUpdate;
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
  controls[8+16+7][3].update = pbDownUpdate;
  controls[8+16+7][3].thresholdPressure = pbThresholdPressure;
  controls[8+16+7][3].maxPressure = pbMaxPressure;
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

struct MenuItem allNotesOffFastMenuItem("notes off", action);
struct MenuItem allNotesOffSlowMenuItem("notes off", allNotesOffSlowAction);
struct MenuItem useUsbMenuItem("USB MIDI", toggle, &useUsbMidi);
struct MenuItem useDinMenuItem("DIN5 MIDI", toggle, &useDinMidi);
struct MenuItem screen10MenuItem("10%", selection, &brightness, 25);
struct MenuItem screen25MenuItem("25%", selection, &brightness, 63);
struct MenuItem screen50MenuItem("50%", selection, &brightness, 127);
struct MenuItem screen75MenuItem("75%", selection, &brightness, 191);
struct MenuItem screen100MenuItem("100%", selection, &brightness, 255);

struct MenuItem mpeHandshakeMenuItem("mpe init", mpeHandshakeAction);

struct MenuItem doVelocityMenuItem("velocity", toggle, &doMpeDynamicVelocity);
struct MenuItem doPressureMenuItem("pressure", toggle, &doMpeDynamicPressure);


struct MenuItem surgeXtPresetMenuItem("Surge XT", selection, &mpeSettings, (uint32_t)&mpeSettingsSurgeXT);
struct MenuItem kspPresetMenuItem("Keystep Pro", selection, &mpeSettings, (uint32_t)&mpeSettingsKSP);
struct MenuItem xv2020PresetMenuItem("XV-2020", selection, &mpeSettings, (uint32_t)&mpeSettingsXV2020);

struct MenuItem outputPresetsMenu("dev presets", submenu, &surgeXtPresetMenuItem, &kspPresetMenuItem, &xv2020PresetMenuItem);
struct MenuItem outputMenu("output", submenu, &useUsbMenuItem, &useDinMenuItem, &mpeHandshakeMenuItem, &outputPresetsMenu);
struct MenuItem controlsMenu("controls", submenu, &doVelocityMenuItem, &doPressureMenuItem);

struct MenuItem screenBrightnessMenu("brightness", submenu, &screen10MenuItem, &screen25MenuItem, &screen50MenuItem, &screen75MenuItem, &screen100MenuItem);
struct MenuItem interfaceMenu("interface", submenu, &screenBrightnessMenu);


struct MenuItem configMenu("setup", submenu, &outputMenu, &controlsMenu, &interfaceMenu);

//struct MenuItem arpMenu("arpeggiator", submenu);

struct MenuItem rootMenu("", submenu, &configMenu, &allNotesOffSlowMenuItem);

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
  ledSetup();
  adcSetup();
  shiftRegisterSetup();
  screenSetup();
  canSetup();
  midiSetup();
  mpeSetup();
  menuSetup();
  statusTextUpdate();
  audioSetup();

  for (int i = 0; i < adcChannels; i++) {
    for (int j = 0; j < adcChannels; j++) {
      if (i==j) {
        calibrationMatrix[i][j] = 100000.0f;
      } else {
        calibrationMatrix[i][j] = 0.0f;
      }
    }
  }

  for (int i = 0; i < maxShiftRegisterBits; i++) {
    for (int j = 0; j < adcChannels; j++) {
      values[i][j] = 0.0;
      resistances[i][j] = 1.0;
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
  //verbose = false;

  if (verbose) {
    uint32_t timestamp = micros();
    uint32_t delta = timestamp - prevTimestamp;

    Serial.print(iteration);
    Serial.print(" ");
    Serial.print(1000.0 / ((float)delta / 1000000.0));
    Serial.print(" midi tx buffer ");
    Serial.println(Serial5.availableForWrite());
    prevTimestamp = timestamp;
  }

  uint32_t prev = 0;
  int32_t focus = 0x7ffffff; /* select a bit to focus on, and don't do any further bit shifting (for debugging) */
  int32_t curBit = -1;

  /* no bit set yet */
  calibrateADCs(false, calibrationMatrix);
  if (verbose) {
    for (int i = 0; i < adcChannels; i++) {
      Serial.print("calibration ");
      for (int j = 0; j < adcChannels; j++) {
        float r = calibrationMatrix[i][j];
        Serial.print(String(1.0/r) + " ");
      }
      Serial.println();
    }
  }

  shiftRegisterClock();
  curBit++;

  while (curBit < maxShiftRegisterBits && curBit <= focus) {
    uint32_t delay = getADCDelay(values[prev], values[curBit]);
    delayMicroseconds(delay);
    int correctionIterations = curBit < 8 ? 4 : 2; /* work harder to get low error on the pot values than the keys */

    if (false || curBit == 5 || curBit == 8) {
      readADCs(verbose, values[curBit], resistances[curBit], calibrationMatrix, correctionIterations);
    } else {
      readADCs(false, values[curBit], resistances[curBit], calibrationMatrix, correctionIterations);
    }

    if (hwversion < 3) {
      for (int j=0; j<adcChannels; j++) {
        values[curBit][j] = deblur(values[prev][j], values[curBit][j]);
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

  prevUsecs = usecs;
  usecs = micros();
  uint32_t delta = usecs - prevUsecs;

  for(int bit = 0; bit < maxShiftRegisterBits; bit++) {
    for (int channel = 0; channel < adcChannels; channel++) {
      auto update = controls[bit][channel].update;
      if (update != nullptr) {
        update(&controls[bit][channel], delta);
      }
    }
  }

  mpeUpdate(delta);

  if (verbose) {
    showResistances();
    //showValues();
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

  if (screenRedrawAge > maxRedrawAge) {
    renderScreen();
    screenRedrawAge = 0;
  }

  if (verbose) {
    Serial.println(String("MIDI Sent: ") + midiMsgsSent + " received: " + midiMsgsReceived);
  }
  iteration++;
}
