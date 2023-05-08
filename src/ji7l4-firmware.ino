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


#define hwversion 1

#include <MIDI.h>

#include <ILI9341_t3.h>
#include <font_Arial.h>
#include <font_ArialBold.h>

#include <ADC.h>
#include <WS2812Serial.h>

/* Pins */

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
#define screenVinPin 3
#define screenDCPin 9

#define touchCSPin 8


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

/* ADCs */

#define adcChannels 4
const int adcPins[adcChannels] = {adc1Pin, adc2Pin, adc3Pin, adc4Pin};
ADC *adc = new ADC();

void adcSetup() {

  pinMode(LED_BUILTIN, OUTPUT);
  for (int i=0; i<adcChannels; i++) {
    int pin = adcPins[i];
    pinMode(pin, INPUT_DISABLE);
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

void readADCs(bool verbose, int* values) {
  /*
  for (int i=0; i < adcChannels; i++) {
    values[i] = adc->adc0->analogRead(adcPins[i]);
  } */

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

  if (verbose) {
    for (int i=0; i<adcChannels; i++) {
      Serial.print(values[i]);
      Serial.print(" ");
    }
    Serial.println("");
  }
}
/* Shift Registers */

void shiftRegisterSetup() {
  pinMode(shiftRegisterOutPin, OUTPUT);
  pinMode(shiftRegisterClockPin, OUTPUT);
}

const int maxShiftRegisterBits = 8+32;

void shiftRegisterReset(int curBit) {
  for (int i=curBit; i<maxShiftRegisterBits; i++) {
    shiftRegisterClock();
  }
  digitalWrite(shiftRegisterOutPin, HIGH);
  delayMicroseconds(1);
  shiftRegisterClock();
  digitalWrite(shiftRegisterOutPin, LOW);
  delayMicroseconds(1);
}

void shiftRegisterClock() {
  digitalWrite(shiftRegisterClockPin, HIGH);
  delayMicroseconds(1);
  digitalWrite(shiftRegisterClockPin, LOW);
  delayMicroseconds(1);
}

/* Screen */

#define TFT_DC      screenDCPin
#define TFT_CS      screenCSPin
#define TFT_RST     255  // 255 = unused, connect to 3.3V
#define TFT_MOSI    sdiPin
#define TFT_SCLK    sckPin
#define TFT_MISO    sdoPin
ILI9341_t3 tft = ILI9341_t3(TFT_CS, TFT_DC, TFT_RST, TFT_MOSI, TFT_SCLK, TFT_MISO);

void screenSetup() {
  pinMode(screenVinPin, OUTPUT);
  analogWrite(screenVinPin, 255);
  tft.begin();
  //tft.setClock(10000000);
  tft.fillScreen(ILI9341_BLUE);
  tft.setTextColor(ILI9341_YELLOW);
  tft.setTextSize(2);
  tft.println("test...");


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

/* MIDI */
MIDI_CREATE_INSTANCE(HardwareSerial, Serial5, dinMidiOut);
//MIDI_CREATE_INSTANCE(HardwareSerial, Serial1, dinMidiIn);
int midiBufferSize = 0;

void midiSetup(){
  dinMidiOut.begin();
  //dinMidiIn.begin();

  midiBufferSize = Serial5.availableForWrite();

  /*
  for (int i = 0; i < 20; i++) {
    dinMidiOut.sendNoteOn(i+30, 100, 1);
    delayMicroseconds(100000);
    dinMidiOut.sendNoteOff(i+30, 100, 1);
    delayMicroseconds(50000);
  } */
}

int midiBufferInUse() {
  return midiBufferSize - Serial5.availableForWrite();
}

bool midiReady() {
  return midiBufferInUse() < 15;
}

bool midiReadyLowPriority() {
  return midiBufferInUse() < 10;
}

#define mpeChannels 4
#define firstMpeChannel 0

#define noOne 0xffff;

struct MpeChannelState{
  int channel;
  bool playing;
  uint32_t age;
  uint8_t lastNote;
  uint8_t lastVolume;
  uint8_t lastFilter;
  uint16_t lastPitchBend;
  uint16_t owner;
  int lastProgramChangeSent;
  uint32_t volumeAge;
};

struct MpeChannelState mpeState[mpeChannels];
uint8_t programChange = 0;

int noteOnCount = 0;
int noteOffCount = 0;

void mpeSetup() {
  for (int channel = firstMpeChannel; channel < mpeChannels; channel++) {
    struct MpeChannelState *state = &mpeState[channel];
    state->channel = channel;
    state->playing = false;
    state->age = 0xffffffff;
    state->lastNote = 0;
    state->lastVolume = 127;
    state->lastFilter = 127;
    state->lastPitchBend = 0;
    state->owner = noOne;
    state->lastProgramChangeSent = -1;
    state->volumeAge = 0xffffffff;
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

struct MpeChannelState *getMpeChannel() {
  int bestChannel = 0;
  uint32_t bestScore = 0;
  if (!midiReady()) {
    return nullptr;
  }
  for (int channel = firstMpeChannel; channel < firstMpeChannel+mpeChannels; channel++) {
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
    dinMidiOut.sendNoteOff(state->lastNote, 127, bestChannel);
    state->playing = false;
    noteOffCount++;
    return state;
  }

  if (state->playing) {
    return nullptr;
  }

  return state;
}

#define middleC 60 /* midi note */
double pbRange = 2.0;

uint8_t pressureCC = 0x01;
//uint8_t pressureCC = 0x07;

double pitchToCents(double pitch) {
  return (log(pitch) / log(2.0)) * 1200.0; 
}

double transpose = 1.0;

struct MpeChannelState *beginMpeNote(double pitch, int owner, double velocity, double pressure) {
  double shift = pitchToCents(pitch * transpose);
  uint8_t volume = pressure * 127.0;
  if (volume > 127) {
    volume = 127;
  }

  if (shift < -(middleC * 100.0)) {
    Serial.print("out of range note (too small) ");
    Serial.println(shift);
    return nullptr;
  } else if (shift > ((128-middleC) * 100.0)) {
    Serial.print("out of range note (too large) ");
    Serial.println(shift);
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
  //uint16_t pb = 0x2000 + ((cents/100.0) * (0x2000/pbRange));

  double range = (double)MIDI_PITCHBEND_MAX - (double)MIDI_PITCHBEND_MIN;
  int pb = (int)((range / (pbRange * 2.0)) * (cents / 100.0));

  if (pb > MIDI_PITCHBEND_MAX) {
    pb = MIDI_PITCHBEND_MAX;
  } else if (pb < MIDI_PITCHBEND_MIN) {
    pb = MIDI_PITCHBEND_MIN;
  }

  state->age = 0;
  state->owner = owner;

  if(midiReady()) {
    state->playing = true;
    state->lastNote = note;
    state->lastPitchBend = pb;
    if (state->lastProgramChangeSent != programChange) {
      dinMidiOut.sendProgramChange(programChange, state->channel+1);
      state->lastProgramChangeSent = programChange;
      Serial.print("sent program change on channel ");
      Serial.println(state->channel);

      //dinMidiOut.sendControlChange(0x0c, 2, state->channel+1);
      Serial.println("set pitchbend range to 2");
    }
    dinMidiOut.sendPitchBend(pb, state->channel+1);

    if (volume != state->lastVolume) {
      dinMidiOut.sendControlChange(pressureCC, volume, state->channel+1);
      state->lastVolume = volume;
      state->volumeAge = 0;
    }

    dinMidiOut.sendNoteOn(note, 127, state->channel+1);
    noteOnCount++;
    Serial.print("sent note-on on channel ");
    Serial.print(state->channel+1);
    Serial.print(" note ");
    Serial.print(note);
    Serial.print(" cents ");
    Serial.println(cents);
  }

  return state;
}

int pressureBackoff = 5000;

void continueMpeNote(struct MpeChannelState *state, double pressure, uint32_t deltaUsecs) {

  /* rate limit pressure updates */
  if (deltaUsecs < pressureBackoff - state->volumeAge) {
    state->volumeAge += deltaUsecs;
    return;
  }

  uint8_t volume = pressure * 127.0;
  if (volume > 127) {
    volume = 127;
  }

  if(midiReadyLowPriority()) {
    dinMidiOut.sendControlChange(pressureCC, volume, state->channel+1);
    state->lastVolume = volume;
    state->volumeAge = 0;
  } else {
    Serial.print(".");
  }
}

bool endMpeNote(struct MpeChannelState *state) {
  if (midiReady()) {
    state->age = 0;
    state->playing = false;
    dinMidiOut.sendNoteOff(state->lastNote, 5, state->channel+1);
    noteOffCount++;
    Serial.print("sent note-off channel ");
    Serial.print(state->channel+1);
    Serial.print(" note ");
    Serial.println(state->lastNote);

    return true;
  }
  return false;
}

/* Setup, Main Loop */

int values[maxShiftRegisterBits][adcChannels];

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

void readPot(int pot, int led){
  if (pot < 1 || pot >> 12) {
    Serial.println("out-of-range pot value");
  }
  int bit, channel;
  if (pot <= 3) {
    bit = 4+pot;
    channel = 2;
  } else {
    bit = pot-3;
    channel = 3;
  }
  int value = values[bit][channel];

  float position = lerp(1900, value, 3750);
  position = pow(position, 4.0);
  

  int color = (clamp8((int)(position*256.0)) << 16);

/*
  if (pot == 4) {
    Serial.print(position);  Serial.print(" "); Serial.println(value);
  }
*/
  leds.setPixel(led, color);
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
  releasing
};

struct Key {
  int index;
  struct Pitch pitch;
  enum KeyState state;
  double intensity;
  struct MpeChannelState *mpeState;
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
  };
  Control(enum ControlType type, int bit, int channel, const char name_[controlNameLen]) : type{type}, bit{bit}, channel{channel} {
    updateFrequency = 1;
    strncpy(name, name_, controlNameLen);
    name[controlNameLen-1] = '\0';
    updateFrequency = 1;
    key = nullptr;
    update = nullptr;
    delay = 0;
  };
  enum ControlType type;
  int bit;
  int channel;
  int updateFrequency;
  uint32_t delay;
  struct Key *key;
  char name[controlNameLen];
  void (*update)(struct Control* control, uint32_t deltaUsecs);
};

#if (hwversion == 0)
int thresholdPressure = 60;
int maxPressure = 400;
#else
int thresholdPressure = 110;
int maxPressure = 600;
#endif

// seconds to fall from full intensity
double releaseRate = 2.0;

void keyUpdate(struct Control* control, uint32_t deltaUsecs) {
  struct Key *key = control->key;
  int value = values[control->bit][control->channel];
  int threshold = thresholdPressure;

  /* debounce */
  if (key->state == playing) {
    threshold -= 3;
  } else {
    threshold += 3;
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

  double intensity = pow(pressure / (double)maxPressure, 0.4);
  double velocity = 1.0;

  double minIntensity = key->intensity - releaseRate * ((double)deltaUsecs/1000000.0);

  if (intensity < minIntensity) {
    intensity = minIntensity;
  }

  if (intensity > 1.0) {
    intensity = 1.0;
  } else if (intensity < 0.0) {
    intensity = 0.0;
  }

  if (intensity > 0.01) {
    switch (key->state) {
      case idle:
      case releasing:
        Serial.print("attempting noteOn ");
        Serial.print(key->pitch.ratio.a);
        Serial.print("/");
        Serial.println(key->pitch.ratio.b);
        key->mpeState = beginMpeNote((double)key->pitch.ratio.a / (double)key->pitch.ratio.b, key->index, velocity, intensity);
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
    }
  } else {
    switch (key->state) {
      case idle:
      case releasing:
        break;
      case playing:
        if (endMpeNote(key->mpeState)) {
          key->mpeState->owner = noOne;
          key->mpeState = nullptr;
          key->state = releasing;
        } else {
          Serial.print(key->index);
          Serial.print(" ");
          Serial.print(intensity);
          Serial.print(" buffer ");
          Serial.print(midiBufferInUse());
          Serial.println(",");
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

bool debounce(struct Control* control, uint32_t deltaUsecs) {
  if (deltaUsecs >= control->delay) {
    control->delay = 0;
    return false;
  }
  
  control->delay -= deltaUsecs;
  return true;
}

#define check_debounce {if (debounce(control, deltaUsecs)) {return;}}

void incProgramChangeUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  if (programChange == 127) {
    return;
  }

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - thresholdPressure;
  if (pressure > 0) {
    programChange++;
    control->delay = 100000;
    Serial.print("programChange ");
    Serial.println(programChange);
  }
}

void decProgramChangeUpdate(struct Control* control, uint32_t deltaUsecs) {
  if (debounce(control, deltaUsecs)) {
    return;
  }

  if (programChange == 0) {
    return;
  }

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - thresholdPressure;
  if (pressure > 0) {
    programChange--;
    control->delay = 100000;
    Serial.print("programChange ");
    Serial.println(programChange);
  }
}

void allNotesOffSlowUpdate(struct Control* control, uint32_t deltaUsecs) {
  if (debounce(control, deltaUsecs)) {
    return;
  }

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - thresholdPressure;
  if (pressure > 0) {
    for (int channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
      for (int note = 0; note < 128; note++) {
        while (!midiReady()) {
          delayMicroseconds(100);
        }
        Serial.println("sending note off");
        dinMidiOut.sendNoteOff(note, 63, channel+1);
      }
    }

    control->delay = 100000;
    Serial.println("all notes off (the slow way)");
  }
}

void allNotesOffUpdate(struct Control* control, uint32_t deltaUsecs) {
  if (debounce(control, deltaUsecs)) {
    return;
  }

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - thresholdPressure;
  if (pressure > 0) {
    for (int channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
      dinMidiOut.sendControlChange(123, 0, channel+1);
    }

    control->delay = 100000;
    Serial.println("all notes off");
  }
}

void transposeUpUpdate(struct Control* control, uint32_t deltaUsecs) {
  if (debounce(control, deltaUsecs)) {
    return;
  }

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - thresholdPressure;
  if (pressure > 0) {
    transpose *= 2.0;
    if (transpose > 8.0) {
      transpose = 8.0;
    }
    control->delay = 200000;
  }
}

void transposeDownUpdate(struct Control* control, uint32_t deltaUsecs) {
  if (debounce(control, deltaUsecs)) {
    return;
  }

  int value = values[control->bit][control->channel];
  int pressure = (4095-value) - thresholdPressure;
  if (pressure > 0) {
    transpose *= 0.5;
    if (transpose < 0.25) {
      transpose = 0.25;
    }
    control->delay = 200000;
  }
}

struct Control controls[maxShiftRegisterBits][adcChannels];

#define CONTROL(type, bit, channel, name) (controls[bit][channel] = Control(type, bit, channel, name))

void controlSetupController() {
  for (int i=0; i<4; i++) {
    for (int j=0; j<adcChannels; j++) {
      controls[i][j] = Control();
    }
  }

  CONTROL(resistor, 0, 0, "cal-10k-r8");
  CONTROL(pressure, 1, 0, "nav-back");
  CONTROL(pressure, 2, 0, "nav-forward");
  CONTROL(pressure, 3, 0, "nav-ok");
  CONTROL(pressure, 4, 0, "nav-cancel");
  CONTROL(pressure, 5, 0, "nav-up");
  CONTROL(pressure, 6, 0, "nav-right");
  CONTROL(pressure, 7, 0, "nav-left");

  CONTROL(resistor, 0, 1, "cal-10k-r13");
  CONTROL(pressure, 1, 1, "nav-down");
  CONTROL(pressure, 2, 1, "button-1");
  CONTROL(pressure, 3, 1, "button-2");
  CONTROL(pressure, 4, 1, "button-4");
  CONTROL(pressure, 5, 1, "button-3");
  CONTROL(pressure, 6, 1, "button-5");
  CONTROL(pressure, 7, 1, "cal-2k-r15");

  CONTROL(resistor, 0, 2, "cal-10k-r17");
  CONTROL(analogIn, 1, 2, "analog-in-1");
  CONTROL(analogIn, 2, 2, "analog-in-2");
  CONTROL(analogIn, 3, 2, "analog-in-3");
  CONTROL(analogIn, 4, 2, "analog-in-4");
  CONTROL(pot, 5, 2, "rv1");
  CONTROL(pot, 6, 2, "rv2");
  CONTROL(pot, 7, 2, "rv3");

  CONTROL(resistor, 0, 3, "cal-10k-r29");
  CONTROL(pot, 1, 3, "rv4");
  CONTROL(pot, 2, 3, "rv5");
  CONTROL(pot, 3, 3, "rv6");
  CONTROL(pot, 4, 3, "rv7");
  CONTROL(pot, 5, 3, "rv8");
  CONTROL(pot, 6, 3, "rv9");
  CONTROL(pot, 7, 3, "rv10");

  controls[2][1].update = incProgramChangeUpdate;
  controls[3][1].update = decProgramChangeUpdate;
  controls[6][1].update = allNotesOffUpdate;
  controls[4][1].update = allNotesOffSlowUpdate;
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

void controlSetupKeybed() {
  int oct = 0;
  CONTROL(resistor, 8, 0, "id0-20k");
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

  CONTROL(pressure, 8+24+1, 0, "4/1");
  CONTROL(pressure, 8+8+6,  3, "control-1");
  CONTROL(pressure, 8+8+7,  3, "control-2");
  CONTROL(pressure, 8+16+6, 3, "spacebar-up");
  CONTROL(pressure, 8+16+7, 3, "spacebar-dn");
  CONTROL(pressure, 8+24+6, 3, "control-4");
  CONTROL(pressure, 8+24+7, 3, "control-3");

  controls[8+8+6][3].update = transposeDownUpdate;
  controls[8+8+7][3].update = transposeUpUpdate;
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

void setup() {
  serialSetup();
  Serial.println("begin setup");
  ledSetup();
  adcSetup();
  shiftRegisterSetup();
  //screenSetup();
  midiSetup();
  mpeSetup();
  controlSetupController();
  controlSetupKeybed();
  delayMicroseconds(100000);
  Serial.println("end setup");
}

uint32_t usecs = 0;
uint32_t prevUsecs = 0;
uint32_t prevTimestamp = 0;

void loop() {
  // put your main code here, to run repeatedly:
  static int iteration = 0;
  bool verbose = (iteration % 4000 == 0);

  if (iteration % 2000 == 0) {
    uint32_t timestamp = micros();
    uint32_t delta = timestamp - prevTimestamp;

    Serial.print(iteration);
    Serial.print(" ");
    Serial.print(2000.0 / ((float)delta / 1000000.0));
    Serial.print(" midi tx buffer ");
    Serial.println(Serial5.availableForWrite());
    prevTimestamp = timestamp;
  }

  int i;
  for (i=0; i<maxShiftRegisterBits; i++) {
    //Serial.print(i);
    //Serial.print(" ");
    shiftRegisterClock();
    delayMicroseconds(10);
    readADCs(false, values[i]);
  }
  shiftRegisterReset(i);

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

  if (verbose)
    showValues();

  if (iteration % 4 == 0) { 
    setLed(1, 7, 0); // nav_left
    setLed(2, 5, 0); // nav_up
    setLed(3, 1, 1); // nav_down
    setLed(4, 6, 0); // nav_right
    //readPot(4, 0);
    //readPot(7, 5);
    int surplus = noteOnCount - noteOffCount;
    if (surplus > 5) {
      surplus = 5;
    }

    leds.setPixel(5, surplus);
    leds.show();
  }

/*
  if (dinMidiIn.read()){
    int note, channel, velocity;
    bool noteCmd = false;

    byte cmd = dinMidiIn.getType();
    switch(cmd) {
      case midi::NoteOn:
        noteCmd = true;
        Serial.print("got NoteOn ");
        break;
      case midi::NoteOff:
        noteCmd = true;
        Serial.print("got NoteOff ");
        break;
      default:
        Serial.println("unexpected midi message");
        break;
    }

    if (noteCmd) {
      note = dinMidiIn.getData1();
      velocity = dinMidiIn.getData2();
      channel = dinMidiIn.getChannel();

      Serial.print("channel "); Serial.print(channel);
      Serial.print(" note "); Serial.print(note);
      Serial.print(" velocity "); Serial.println(velocity);
    }
  }
  */
  iteration++;
}
