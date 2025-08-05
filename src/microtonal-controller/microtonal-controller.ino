/**@file microtonal-controller.ino */

/*
Copyright 2023-2025 Jim Snow, Desiderata Systems LLC

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

#define fwversion "1.1.1"

#define hwversion 4

#include <MIDI.h>

#include <ILI9341_t3.h>
#include <font_Arial.h>
#include <font_ArialBold.h>

#include <ADC.h>
#include <WS2812Serial.h>
#include <FlexCAN_T4.h>

/* Ownership */
#define noOne 0xffff

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

#if (hwversion >= 3)
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

/* hwversion 4 has no software-visible changes from hwversion 3 */

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

/* ADCs */

#define adcChannels 4

float zeroPressureResistance = 5500.0f;
float maxPressureResistance = 800.0f;

/*
 * Force can be outside range of 0 (minimum force) to 1 (maximum force)
 */
 
float resistanceToForce(float r, float area = 1.0f) {
  return 1.0f - lerpNoClamp(max(area, 2.0f) / maxPressureResistance, 1.0f / r, area / zeroPressureResistance);
}

const int maxShiftRegisterBits = 8+32;

/* Screen */


/* Menu */

bool lock = false;  /* disables most controls when true */

enum menuItemType {
  action,
  toggle,
  value,
  floatValue,
  selection,
  submenu,
  empty
};

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

void statusTextUpdate();
void status1TextUpdate(String text, uint32_t usecs = 0);

#include "menu.h"

uint32_t brightness = 127;
uint32_t brightnessSet = 255;

struct MenuItem emptyMenuItem = MenuItem("", empty);

/* CAN BUS */


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
uint32_t pressureBackoff = 5000;

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

/* used for doing "all notes off" the manual way */
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

// monitor serial send buffer to avoid overrunning it
// (we should probably do likewise for usb)
int midiBufferInUse() {
  if (useDinMidi) {
    int avail = Serial5.availableForWrite();
    int inUse = midiBufferSize - avail;
    if (inUse < 0) {
      return 0;
    }
    return inUse;
  }
  return 0;
}

// test if there's enough send buffer space to send at least a few more MIDI messages
bool midiReady() {
  return midiBufferInUse() < 15;
}

// test if the send buffer backlog is short enough we can send low-priority messages
bool midiReadyLowPriority() {
  return midiBufferInUse() < 10;
}

// test if the send buff is empty
bool midiIdle() {
  return midiBufferInUse() == 0;
}

// wait for send buffer space to become available
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

/* MPE */
float mpePolyAfterTouchMin = 0.75f;
float mpePolyAfterTouchMax = 1.25f;

uint32_t mpeChannels = 16;
uint32_t mpeChannelsSet = mpeChannels;
uint32_t firstMpeChannel = 0;

struct MpeChannelState{
  int channel;
  bool playing;
  uint32_t age;
  uint8_t lastNote;
  uint8_t lastVolume;
  uint8_t lastPolyAT;
  float volume;
  float filter;
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

bool doCCPassThrough = true;

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

  score += state->age >> 9;

  score += state->lastVolume << 23;

  return score;
}

bool skipChannel10 = false;

struct MpeChannelState *getMpeChannel() {
  uint32_t bestChannel = 0;
  uint32_t bestScore = 0;
  if (!midiReady()) {
    return nullptr;
  }
  for (uint32_t channel = firstMpeChannel; channel < firstMpeChannel+mpeChannels; channel++) {
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
    endMpeNote(state, state->owner);
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
uint32_t minMpePressure = 0;

double pbUp = 0.0;
double pbDown = 0.0;

/* seconds from full to none, or vice versa */
double attack = 0.0;
double decay = 1.0;

// 1/seconds to fall from full intensity
float volumeReleaseRate = 2.0f;
float volumeAttackRate = 1000000.0f;
float filterReleaseRate = 2.0f;
float filterAttackRate = 1000000.0f;

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
bool doMpePressureFilter = false;
float mpeStaticVelocity = 0.75;
//float minVelocity = 0.2;
bool doMpeChannelPressure = false;

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
double masterPbDownRange = 1.0/2.0;  // 2.0/3.0;

uint32_t masterPbAge = 0;
int16_t lastMasterPb = 0;
uint32_t masterPbBackoff = 8000;

double masterPbRange = 12.0;

bool forcePbRange = false;

bool delayNoteOff = false;
bool noteOnFirst = false;
bool doFB01Setup = false;

uint32_t minVelocity = 1; /* minimum MIDI note-on velocity, range: 1-127 */

double pitchReferenceHz() {
  static double c = 440.0 / (pow(pow(2, 1.0/12.0), 9));

  return c * transpose;
}

double calculatePitchBendCents(double pbUp, double pbDown) {
  double pb = clampDoubleBipolar(pbUp - pbDown);
  double cents = pb > 0
    ? pitchToCents(masterPbUpRange) * pb
    : pitchToCents(masterPbDownRange) * -pb;

  return cents;
}

double calculatePitchBend(double pbUp, double pbDown) {
  return centsToPitch(calculatePitchBendCents(pbUp, pbDown));
}

int16_t calculateMidiPitchBend(double pbUp, double pbDown, double pitch, double pbRange) {

  double cents = calculatePitchBendCents(pbUp, pbDown) + pitchToCents(pitch);
    
  double range = (double)MIDI_PITCHBEND_MAX - (double)MIDI_PITCHBEND_MIN;
  int pbi = (int)((range / (pbRange * 2.0)) * (cents / 100.0));
  if (pbi > MIDI_PITCHBEND_MAX) {
    return MIDI_PITCHBEND_MAX;
  } else if (pbi < MIDI_PITCHBEND_MIN) {
    return MIDI_PITCHBEND_MIN;
  }

  return pbi;
}

void fb01SysEx1Byte(uint8_t channel, uint8_t param, uint8_t data) {
  uint8_t system = 0;
  uint8_t instrument = channel;

  uint8_t buffer[8];
  buffer[0] = 0xf0; // status
  buffer[1] = 0x43; // Yamaha
  buffer[2] = 0x75; // sub-status
  buffer[3] = system & 0x0f; // which system
  buffer[4] = instrument | 0x18; // which instrument
  buffer[5] = param & 0x7f; // parameter
  buffer[6] = data & 0x7f;  // 7 bit payload
  buffer[7] = 0xf7; // end

  midiReadyWait();
  midiSysEx(8, buffer);
}

void fb01SysEx2Byte(uint8_t channel, uint8_t param, uint8_t data) {
  uint8_t system = 0;
  uint8_t instrument = channel;

  uint8_t buffer[9];
  buffer[0] = 0xf0; // status
  buffer[1] = 0x43; // Yamaha
  buffer[2] = 0x75; // sub-status
  buffer[3] = system & 0x0f;
  buffer[4] = instrument | 0x18;
  buffer[5] = param & 0x7f;
  buffer[6] = data & 0x0f;
  buffer[7] = (data >> 4) &0x0f; 
  buffer[8] = 0xf7;

  midiReadyWait();
  midiSysEx(9, buffer);
}

void setPitchBendRange(int channel) {
  midiReadyWait();

  if (doFB01Setup) {
    fb01SysEx1Byte(channel, 0x0c, pbRange); // set pitch bend
    fb01SysEx1Byte(channel, 0x0a, 0);       // turn off lfo
    fb01SysEx1Byte(channel, 0x06, 0);       // turn off detune
    fb01SysEx2Byte(channel, 0x7b, 0x20 | (pbRange & 0x0f));  // set pitch bend control to pitch wheel, and set bend range 
  } else {
    midiRPN(0, 0, pbRange, channel + 1);
  }
}

void prepareChannel(struct MpeChannelState *state) {
  int channel = state->channel;

  /*
   * If we're updating the LSB, we have to send the MSB first
   * even if it hasn't changed. (On the XV-2020 at least.)
   */
  if (doMsb && (state->lastBankMsb != mpeBankMsb || state->lastBankLsb != mpeBankLsb)) {
    midiReadyWait();
    midiControlChange(0, mpeBankMsb, channel + 1);
    state->lastBankMsb = mpeBankMsb;
    state->lastProgramChangeSent = 128;  /* force new program change to be sent */
    Serial.println("sent msb " + String(mpeBankMsb) + " channel " + String(channel));
  }

  if (doLsb && state->lastBankLsb != mpeBankLsb) {
    midiReadyWait();
    if (doFB01Setup) {
      fb01SysEx1Byte(state->channel, 0x04, mpeBankLsb);
    } else {
      midiControlChange(32, mpeBankLsb, channel + 1);
    }
    state->lastBankLsb = mpeBankLsb;
    state->lastProgramChangeSent = 128; /* force new program change to be sent */
    Serial.println("sent lsb" + String(mpeBankLsb) + " channel " + String(channel));
  }

  if (doProgramChange && state->lastProgramChangeSent != programChange) {
    midiReadyWait();
    midiProgramChange(programChange, channel + 1);
    state->lastProgramChangeSent = programChange;

    // set pitch bend range
    if (forcePbRange) {
      setPitchBendRange(state->channel);
    }
    Serial.println("sent program change " + String(programChange) + " channel " + String(channel));
  }
}

bool bendUpOnly = false;
bool bendDownOnly = false;
bool enableBender = true;

struct MpeChannelState *beginMpeNote(double pitch, double velocity, double pressure, uint16_t owner) {

  uint32_t v = (int) (1.0 + (velocity*126));
  if (v > 127) {
    v = 127;
  } else if (v < minVelocity) {
    v = minVelocity;
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

  if (bendDownOnly && cents > 0.0f) {
    semitones++;
    cents -= 100.0f;
  }

  uint32_t note = middleC + semitones;
  int16_t pb;
  if (midiType == mpe) {
    pb = calculateMidiPitchBend(0.0, 0.0, centsToPitch(cents), (double)pbRange);
  } else {
    pb = calculateMidiPitchBend(pbUp, pbDown, centsToPitch(cents), (double)pbRange);
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

  if(midiReady()) {
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
      state->playing = true;
    }

    continueMpeNote(state, pressure, 0, owner);

    if (!noteOnFirst) {
      midiNoteOn(note, v, state->channel+1);
      state->playing = true;
    }

    Serial.println("note-on channel" + String(state->channel+1) + " note " + String(note) + " cents " + String(cents) + " (" + String(centsToPitch(cents)) + ") pb " + String(pb) + " pressure " + pressure);
  }

  state->age = 0;

  return state;
}

void continueMpeNote(struct MpeChannelState *state, double pressure, uint32_t deltaUsecs, uint16_t owner) {

  if (state == nullptr || state->owner != owner) {
    return;
  }

  uint32_t concurrentNotes = noteOnCount - noteOffCount;

  float volume, filter;

  if (pressure < 0.0f) {
    pressure = 0.0f;
  }

  float minVolume = state->volume - (volumeReleaseRate * (state->volume + 0.2f)) * ((float)deltaUsecs/1000000.0f);
  float maxVolume = state->volume + volumeAttackRate  * ((float)deltaUsecs/1000000.0f) * (pressure);
  if (pressure < minVolume) {
    volume = minVolume;
  } else if (pressure > maxVolume) {
    volume = maxVolume;
  } else {
    volume = pressure;
  }
  state->volume = volume;

  float minFilter = state->filter - (filterReleaseRate * (state->volume + 0.2f)) * ((float)deltaUsecs/1000000.0f);
  float maxFilter = state->filter + filterAttackRate  * ((float)deltaUsecs/1000000.0f) * (pressure * pressure);
  if (pressure < minFilter) {
    filter = minFilter;
  } else if (pressure > maxFilter) {
    filter = maxFilter;
  } else {
    filter = pressure;
  }
  state->filter = filter;

  /* rate limit pressure updates */
  if (state->volumeAge + deltaUsecs < pressureBackoff * concurrentNotes) {
    state->volumeAge += deltaUsecs;
    return;
  }

  int max = maxMpePressure;
  uint16_t vol = max;

  if (doMpeDynamicPressure) {
    vol = minMpePressure + (volume * (float)(maxMpePressure - minMpePressure));
    if (vol > max) {
      vol = max;
    }
  }

  if(midiReadyLowPriority()) {
    if (volume != state->lastVolume) {
      if (doMpeChannelPressure) {
        midiAfterTouch(vol, state->channel+1);
        //Serial.println("aftertouch " + String(volume));
      } else {
        if (pressureCC != 128) {
          midiControlChange(pressureCC, vol, state->channel+1);
        }
      }
      state->lastVolume = vol;
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

    if (doMpePressureFilter) {
      uint8_t fc = (uint16_t)(filter * 127.0f);
      if (fc != state->lastFilter) {
        midiControlChange(74, fc, state->channel+1);
        state->lastFilter = fc;
      }

    }
  } else {
    Serial.print(".");
    state->volumeAge += deltaUsecs;
  }

  if (state->playing && delayNoteOff && volume <= 0.0f) {
    endMpeNote(state, owner);
  }
  
  return;
}

void endMpeNote(struct MpeChannelState *state, uint16_t owner) {
  if (state == nullptr || state->owner != owner) {
    return;
  }

  state->age = 0;
  state->playing = false;
  midiReadyWait();
  midiNoteOff(state->lastNote, 63, state->channel+1);

  Serial.println("sent note-off channel " + String(state->channel+1) + " note " + String(state->lastNote));
  state->lastPolyAT = 0;
  state->owner = noOne;
}

void doMpeMasterPitchbend(double pbUp, double pbDown, uint32_t deltaUsecs) {
  if (!enableBender) {
    pbUp = 0.0; pbDown = 0.0;
  }

  if (deltaUsecs < masterPbBackoff - masterPbAge) {
    masterPbAge += deltaUsecs;
    return;
  }

  if (midiType == mpe || midiType == tuningtable) {
    int16_t pbi = calculateMidiPitchBend(pbUp, pbDown, 1.0, masterPbRange);
    if (pbi != lastMasterPb) {
      midiPitchBend(pbi, 1);
      lastMasterPb = pbi;
    }
  } else if (midiType == multitimbral || midiType == monotimbral || midiType == monovoice) {
    for (uint32_t channel = firstMpeChannel; channel < firstMpeChannel+mpeChannels; channel++) {
      if (mpeState[channel].playing) {
        int16_t pbi = calculateMidiPitchBend(pbUp, pbDown, mpeState[channel].lastBendInterval, (double)pbRange);
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


struct MidiTuningTableEntry* beginTuningTableNote(double pitch, double velocity, double pressure, uint16_t owner) {
  if (midiTuningTable == nullptr || midiTuningTableSize == 0) {
    return nullptr;
  }

  struct MidiTuningTableEntry *tte = tuningTableLookup(pitch, midiTuningTable, midiTuningTableSize);
  
  int v = (int) (1.0 + (velocity*126));
  if (v > 127) {
    v = 127;
  } else if (v < 1) {
    v = 1;
  }

  midiReadyWait();
  midiNoteOn(tte->noteNumber, v, tte->channel + 1);
  tte->users++;
  tte->lastPressure = 127;
  tte->pressureAge = 0;
  return tte;
}

void continueTuningTableNote(struct MidiTuningTableEntry *tte, double pressure, uint32_t deltaUsecs, int16_t owner) {
  uint32_t concurrentNotes = noteOnCount - noteOffCount;

  if (pressure < 0.0f) {
    pressure = 0.0f;
  }

  /* rate limit pressure updates */
  if (tte->pressureAge + deltaUsecs < pressureBackoff * concurrentNotes) {
    tte->pressureAge += deltaUsecs;
    return;
  }

  if (!midiReady()) {
    return;
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
}

void endTuningTableNote(struct MidiTuningTableEntry* tte, uint16_t owner) {
  if (midiReady()) {
    midiNoteOff(tte->noteNumber, 63, tte->channel + 1);
    tte->users--;
  }
}

/* return the buffer offset for a given midi note and sub-byte */
#define eTuningTableOffset(index, field) (7 + (index * 4) + field)

/*
 * Send tuning table for the Grey Matter E! expansion board for Yamaha DX7
 * (midiTuningTable should already be initialized before calling this)
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

bool doLocalSynth = true;

struct VoiceHandle {
  enum midiType midiType;
  union {
    struct MpeChannelState *mpeState;
    struct MidiTuningTableEntry *tte;
  };
  struct SynthVoice *voice;
};

void beginNote(struct VoiceHandle &voiceHandle, float pitch, float velocity, float pressure, uint16_t owner) {
  switch (midiType) {
    case tuningtable:
      {
        struct MidiTuningTableEntry *tte = beginTuningTableNote(pitch, velocity, pressure, owner);
        voiceHandle.tte = tte;
        break;
      }
    default:
      {
        struct MpeChannelState *state = beginMpeNote(pitch, velocity, pressure, owner);
        voiceHandle.mpeState = state;
        break;
      }
  }

  if (doLocalSynth) {
    voiceHandle.voice = beginSynthNote(pitch, velocity, pressure, owner);
  }

  voiceHandle.midiType = midiType;
}

void continueNote(struct VoiceHandle& voiceHandle, float pressure, uint32_t deltaUsecs, uint16_t owner) {
  switch (voiceHandle.midiType) {
    case tuningtable:
      continueTuningTableNote(voiceHandle.tte, pressure, deltaUsecs, owner);
    default:
      continueMpeNote(voiceHandle.mpeState, pressure, deltaUsecs, owner);
  }

  if (doLocalSynth) {
    continueSynthNote(voiceHandle.voice, pressure, deltaUsecs, owner);
  }
}

void endNote(struct VoiceHandle& voiceHandle, uint16_t owner) {
  if (doLocalSynth) {
    endSynthNote(voiceHandle.voice, owner);
  }
  switch (voiceHandle.midiType) {
    case tuningtable:
      endTuningTableNote(voiceHandle.tte, owner);
    default:
      endMpeNote(voiceHandle.mpeState, owner);
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
#define multicastTimbreFlag (1 << 10) /* CCs sent to mpe multicast channel don't propagate to others, so do it manually */
#define delayNoteOffFlag    (1 << 11) /* defer note-off according to pressure slew limiter */
#define noteOnFirstFlag     (1 << 12) /* send note on before pressure */
#define maxPressure126Flag  (1 << 13) /* limit max pressure to 126 */
#define bendUpOnlyFlag      (1 << 14) /* don't bend down for pitch correction, only up */
#define bendDownOnlyFlag    (1 << 15) /* don't bend up for pitch correction, only down */
#define useETuningTableFlag (1 << 16) /* use Grey Matter E! tuning table format */ 
#define doFB01SetupFlag     (1 << 17) /* FB01 has some specific sysex configuration */
#define minVelocity16Flag   (1 << 18) /* set min velocity to 16 */
#define minVelocity32Flag   (1 << 19) /* set min velocity to 32 */
#define minPressure32Flag   (1 << 20) /* set min pressure to 32 */
#define pressureFilterFlag  (1 << 21) /* modulate the filter in addition to volume with key pressure */
#define channelPressureFlag (1 << 22) /* use channel pressure rather than CC for key pressure (mpe default) */

struct MpeBank {
  uint8_t lsb;
  uint8_t msb;
  String desc;
};

struct MpeBankMap {
  uint8_t bankLength;
  struct MpeBank* banks;
};

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
  struct MpeBankMap *bankMap;
  uint32_t flags;
};

const struct MpeSettings mpeSettingsUsbMidi = {multitimbral, 16,  2,  1000, 7,   0,  127, 0,  0,  127, 0,  nullptr, useUsbFlag | polyAfterTouchFlag};
const struct MpeSettings mpeSettingsXV2020  = {multitimbral, 16,  2, 15000, 7,   87, 87,  87, 64, 67,  64, nullptr, useDinFlag | skipChannel10Flag | gmFlag | gm2Flag | noPressureFlag};
const struct MpeSettings mpeSettingsRD300NX = {multitimbral, 16,  2,  1000, 7,   87, 121, 87, 0,  127, 0,  nullptr, useDinFlag | gmFlag};
const struct MpeSettings mpeSettingsFB01    = {multitimbral, 8,   2,  1000, 7,   0,  0,   0,  0,  6,   0,  nullptr, useDinFlag | doFB01SetupFlag | forcePbRangeFlag};
const struct MpeSettings mpeSettingsKSP     = {multitimbral, 4,   2,  1000, 1,   0,  0,   0,  0,  0,   0,  nullptr, useDinFlag};
const struct MpeSettings mpeSettingsMT2     = {multitimbral, 4,   2,  1000, 128, 0,  0,   0,  0,  0,   0,  nullptr, useDinFlag | pressureFilterFlag | delayNoteOffFlag | channelPressureFlag};
const struct MpeSettings mpeSettingsTrinity = {multitimbral, 16,  2,  1000, 7,   0,  0,   0,  0,  3,   0,  nullptr, useDinFlag | forcePbRangeFlag | gmFlag};
const struct MpeSettings mpeSettingsX50     = {multitimbral, 16,  2,  1000, 7,   63, 63,  63, 0,  3,   0,  nullptr, useDinFlag | forcePbRangeFlag | gmFlag};
const struct MpeSettings mpeSettingsSP300   = {multitimbral, 16,  2,  1000, 7,   0,  0,   0,  0,  0,   0,  nullptr, useDinFlag | skipChannel10Flag | minVelocity16Flag};
const struct MpeSettings mpeSettingsSurgeXT = {mpe,          16, 48,  1000, 128, 0,  0,   0,  0,  0,   0,  nullptr, useUsbFlag | multicastTimbreFlag | maxPressure126Flag };
const struct MpeSettings mpeSettingsProteus = {multitimbral, 16,  2,  1000, 7,   0,  4,   4,  0,  7,   0,  nullptr, useDinFlag | gmFlag | forcePbRangeFlag | minVelocity32Flag | minPressure32Flag};
const struct MpeSettings mpeSettingsProteusDemo = {multitimbral, 16,  2,  1000, 7,   0,  4,  80,  0,  7,   0,  nullptr, useDinFlag | gmFlag | forcePbRangeFlag | minVelocity32Flag | minPressure32Flag};
const struct MpeSettings mpeSettingsMox8    = {multitimbral, 16,  2,  1400, 7,   63, 63,  63, 0,  7,   0,  nullptr, useDinFlag | gmFlag | forcePbRangeFlag | ccResetFlag};
const struct MpeSettings mpeSettingsPhatty  = {monovoice,    1,   2,  1000, 19,  0,  0,   0,  0,  0,   0,  nullptr, useDinFlag | forcePbRangeFlag};
const struct MpeSettings mpeSettingsDx7E    = {tuningtable,  1,   1,  1000, 128, 0,  0,   0,  0,  0,   0,  nullptr, useDinFlag | useETuningTableFlag | noPressureFlag};
const struct MpeSettings mpeSettingsDexed   = {tuningtable,  1,   1,  1000, 128, 0,  0,   0,  0,  0,   0,  nullptr, useUsbFlag | noPressureFlag };
const struct MpeSettings mpeSettingsQSynth  = {multitimbral, 16,  2,  1000, 11,  0,  0,   0,  0,  0,   0,  nullptr, useUsbFlag | skipChannel10Flag | minVelocity32Flag};


void sendMpeZones() {
  midiReadyWait();
  midiRPN(0, 6, mpeChannels, 1);
}

const struct MpeSettings *currentMpeSettings = nullptr;
const struct MpeSettings *mpeSettings = nullptr;

void applyMpeSettings(const struct MpeSettings *settings) {
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
  noteOnFirst = (settings->flags & noteOnFirstFlag) != 0;
  bendUpOnly = (settings->flags & bendUpOnlyFlag) != 0;
  bendDownOnly = (settings->flags & bendDownOnlyFlag) != 0;
  doFB01Setup = (settings->flags & doFB01SetupFlag) != 0;
  minVelocity = (settings->flags & minVelocity16Flag) != 0
    ? 16
    : (settings->flags & minVelocity32Flag) != 0 ? 32 : 1;
  minMpePressure = (settings->flags & minPressure32Flag) != 0 ? 32 : 0;
  doMpePressureFilter = (settings->flags & pressureFilterFlag) != 0;
  doMpeChannelPressure = (settings->flags & channelPressureFlag) != 0 || settings->midiType == mpe;


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

  if (doFB01Setup) {
    for (int i = 0; i < 8; i++) {
      fb01SysEx1Byte(i, 0x01, i);       // set instrument i to midi channel i
      fb01SysEx1Byte(i, 0x00, 1);       // 1 voice per instrument
      fb01SysEx1Byte(i, 0x06, 0);       // detune
      fb01SysEx1Byte(i, 0x07, 2);       // octave transpose (range 0-4, with 2 as center)
      fb01SysEx1Byte(i, 0x09, 64);      // pan middle
      fb01SysEx1Byte(i, 0x0a, 0);       // disable LFo
      fb01SysEx1Byte(i, 0x0c, pbRange); // pitchbend range
    }
  }

  currentMpeSettings = settings;
  mpeSettings = settings;

  statusTextUpdate();
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
      endMpeNote(state, state->owner);

      midiReadyWait();
      if (doMpeChannelPressure) {
        midiAfterTouch(1.0f, i+1);
      } else {
        midiControlChange(pressureCC, 1.0f, i+1);
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
  for (uint32_t channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
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
  }

  applyMpeSettings(&mpeSettingsSurgeXT);
  //applyMpeSettings(&mpeSettingsMT2);
  //applyMpeSettings(&mpeSettingsProteusDemo);
  //programChange = 17;
}

/* Used by MOX8 device preset */
void mpeCCReset(){
  for (uint32_t channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
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

uint32_t idleCount = 0;

void mpeUpdate(uint32_t deltaUsecs) {
  uint16_t channelMask = 0;
  for (uint32_t channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
    struct MpeChannelState *state = &mpeState[channel];
    state->age += deltaUsecs;
    if (state->owner == noOne) {
      if (state->playing) {
        continueMpeNote(state, 0.0f, deltaUsecs, noOne);
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
        for (uint32_t channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels && !congested; channel++) {
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

float reverseLerp(float a, float b, float c) {
  float range = c-a;
  float value = a + b*range;
  if (value < a) {
    return a;
  } else if (value > c) {
    return c;
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
  ratioType,
  frequencyType,
  centsType,
  edoType
};

struct RatioPitch {
  int a;
  int b;
};

struct EdoPitch {
  int edoBase;
  int steps;
};

int gcd_(int a, int b) {
  if (a == b) {
    return a;
  } else if (a < b) {
    return gcd_(b, a);
  } else {
    return gcd_(a-b, b);
  }
}

int gcd(int a, int b) {
  Serial.print("gcd(" + String(a) + "," + String(b) + ") = ");
  int ret = gcd_(a,b);
  Serial.println(String(ret));
  return ret;
}

struct Pitch {
  void simplify() {
    if (type != ratioType) {
      return;
    }

    int div = gcd(ratio.a, ratio.b);

    ratio.a /= div;
    ratio.b /= div;
  }
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
  struct Pitch originalPitch;
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
  float max = 8500.0;
  float min = 600.0;

  float valueUpperBound = 100.0;
  float valueLowerBound = 0.0;
  float hysteresis = 800.0;
  uint32_t data = 0;
  float initial = 0.0f;
  enum knobState state = disabled;
  bool activeByDefault = false;
  void (*action)(uint32_t data, float current) = nullptr;
  String label = "";
};

float filterAmount = 0.5f;
float resonanceAmount = 0.1;
float pressureExponent = 0.8f;

void knobVolumeAction(uint32_t cc, float value) {
  setSynthVolume(value);
  knobCCAction(cc, value);
}

void knobModAction(uint32_t cc, float value) {
  setSynthModulation(value);
  knobCCAction(cc, value);
}

void knobFilterAction(uint32_t cc, float value) {
  filterAmount = value;
  knobCCAction(cc, value);
}

void knobResonanceAction(uint32_t cc, float value) {
  resonanceAmount = value;
  knobCCAction(cc, value);
}

void knobReverbAction(uint32_t cc, float value) {
  setSynthReverb(value);
  knobCCAction(cc, value);
}

void knobCCAction(uint32_t cc, float value) {
  if (cc == 74 && doMpePressureFilter) {
    return;
  }
  mpeMulticastCC(cc, value);
}

void knobVolumeReleaseRateAction(uint32_t unused, float value) {
  value = audioTaper(value);
  if (value < 0.0001f) {
    value = 0.0001f;
  }

  volumeReleaseRate = 0.3f / value;
}

void knobVolumeAttackRateAction(uint32_t unused, float value) {
  value = audioTaper(value);

  if (value < 0.0001f) {
    value = 0.0001f;
  }

  volumeAttackRate = 1.0f / value;
}

void knobFilterReleaseRateAction(uint32_t unused, float value) {
  value = audioTaper(value);
  if (value < 0.0001f) {
    value = 0.0001f;
  }

  filterReleaseRate = 0.3f / value;
}

void knobFilterAttackRateAction(uint32_t unused, float value) {
  value = audioTaper(value);
  if (value < 0.0001f) {
    value = 0.0001f;
  }

  filterAttackRate = 1.0f / value;
}

void knobPressureExponentAction(uint32_t unused, float value) {
  value = clamp(value);

  float min = 0.25f;
  float max = 1.5f;

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

void keyUpdate(struct Control* control, uint32_t deltaUsecs) {
  struct Key *key = control->key;
  float pressure = forces[control->bit][control->channel];

  float lastPressure = key->lastPressure;
  key->lastPressure = pressure;

  float velocity = mpeStaticVelocity;

  if (doMpeDynamicVelocity) {
    float delta = pressure - lastPressure;
    velocity = delta > 0
      ? pow( ( delta * 8000.0f) / (float)deltaUsecs, 0.3f)
      : pow( (-delta * 8000.0f) / (float)deltaUsecs, 0.3f);
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

  key->intensity = intensity >= 0.0f ? intensity : 0.0f;

  if (pressure > 0.01f && (lastPressure > 0.01f || pressure >= 1.0f || !doMpeDynamicVelocity)) {
    switch (key->state) {
      case idle:
      case releasing:
        beginNote(key->voiceHandle, (double)key->pitch.ratio.a / (double)key->pitch.ratio.b, velocity, intensity, key->index);
        Serial.println("noteOn " + String(key->pitch.ratio.a) + "/" + String(key->pitch.ratio.b) + " v=" + String(velocity));
        key->state = playing;
        break;
      case playing:
        continueNote(key->voiceHandle, intensity, deltaUsecs, key->index);
        break;
    }
  } else if (pressure <= 0.0f && (lastPressure <= 0.0f || !doMpeDynamicVelocity)) {
    switch (key->state) {
      case idle:
      case releasing:
        break;
      case playing:
        if (delayNoteOff && doMpeDynamicPressure) {
          continueNote(key->voiceHandle, intensity, deltaUsecs, key->index);
        } else {
          continueNote(key->voiceHandle, intensity, deltaUsecs, key->index);
          endNote(key->voiceHandle, key->index);
          key->state = releasing;
        }
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
  if (lock) {
    return;
  }

  int knobIndex = control->data;
  struct Knob *knob = &knobs[knobIndex];
  float r = resistances[control->bit][control->channel];

  if (r > knob->max) {
    r = knob->max;
  }

  if (r < knob->min) {
    r = knob->min;
  }

  float h = (knob->hysteresis * 0.5f) + (knob->hysteresis * (r / knob->max) * 1.0f);  /* readings are less stable in the upper range */

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

  bool ascending = false;
  bool descending = false;

  if (r > knob->valueUpperBound) {
    knob->valueUpperBound = max(knob->valueUpperBound, r + (h * 0.05f));
    knob->valueLowerBound = max(knob->valueLowerBound, r - (h * 0.95f));
    ascending = true;
  } else if (r < knob->valueLowerBound) {
    knob->valueLowerBound = min(knob->valueLowerBound, r - (h * 0.05f));
    knob->valueUpperBound = min(knob->valueUpperBound, r + (h * 0.95f));
    descending = true;
  }

  if (ascending || descending) {
    if (ascending) {
      r = (-h * 0.5) + r;
    } else {
      r = (h * 0.5) + r;
    }

    float value = clamp((r - knob->min) / (knob->max - knob->min));

    Serial.println("knob " + String(knobIndex + 1) + " value " + String(value));

    if (knob->action != nullptr) {
      knob->action(knob->data, value);
    }

    status1TextUpdate(knob->label, 500000);
  }
}

void scrollUpUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  float force = forces[control->bit][control->channel];
  if (force > 0.0f) {
    menuScroll(-1);
    control->delay = (uint32_t)((1.2f - clamp(force)) * 200000.0f);
  }
}

void scrollDownUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  float force = forces[control->bit][control->channel];
  if (force > 0.0f) {
    menuScroll(1);
    control->delay = (uint32_t)((1.2f - clamp(force)) * 200000.0f);
  }
}

void allNotesOffFast() {
  for (uint32_t channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
    while (!midiReady()) {
      delayMicroseconds(100);
    }
    midiControlChange(123, 0, channel+1);
  }
}

void allNotesOffSlow() {
  for (uint32_t channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
    for (uint32_t note = 0; note < 128; note++) {
      while (!midiReady()) {
         delayMicroseconds(100);
      }
      midiNoteOffNoCount(note, 63, channel+1);

      /* even on USB it seems we need some delay in here */
      delayMicroseconds(100);
    }
  }

  Serial.println("all notes off");
}

void resetAllControllers() {
  for (uint32_t channel = firstMpeChannel; channel < firstMpeChannel + mpeChannels; channel++) {
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

  if (lock) {
    output = output + " L";
  }

  setStatusText(" " + name + String(octave+4) + " " + type + output);
}

void transposeUpdate(struct Control* control, uint32_t deltaUsecs, float transposeAmount) {
  check_debounce;

  float maxTranspose = lock ? 2.0f : 8.0f;
  float minTranspose = lock ? 0.5f : 0.25f;

  float force = forces[control->bit][control->channel];
  if (force > 0.02f) {
    if (!control->held) {
      transpose *= transposeAmount;
      if (transpose > maxTranspose) {
        transpose = maxTranspose;
      } else if (transpose < minTranspose) {
        transpose = minTranspose;
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
  if (lock) {
    return;
  }
  transposeUpdate(control, deltaUsecs, semitone);
}

void transposeDownSemitoneUpdate(struct Control* control, uint32_t deltaUsecs) {
  if (lock) {
    return;
  }
  transposeUpdate(control, deltaUsecs, 1.0f/semitone);
}

void editValueIncrementUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  float force = forces[control->bit][control->channel];
  if (force > 0.0f) {
    incrementMenuValue();
    control->delay = (uint32_t)((1.0 - force) * 200000.0);
  }
}

void editValueDecrementUpdate(struct Control* control, uint32_t deltaUsecs) {
  check_debounce;

  float force = forces[control->bit][control->channel];
  if (force > 0.0f) {
    decrementMenuValue();
    control->delay = (uint32_t)((1.0 - force) * 200000.0);
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
  controls[2][0].update = menuButtonUpdate;
  controls[2][0].data = fwdText;
  controls[3][0].update = menuButtonUpdate;
  controls[3][0].data = okText;
  controls[4][0].update = menuButtonUpdate;
  controls[4][0].data = cancelText;

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
  knobs[3].action = &knobModAction;
  knobs[3].state = active;
  knobs[3].label = "mod";

  knobs[7].data = 74; /* filter cutoff / MPE timbre */
  knobs[7].action = &knobFilterAction;
  knobs[7].state = active;
  knobs[7].label = "timbre";

  knobs[8].data = 71; /* filter resonance */
  knobs[8].action = &knobResonanceAction;
  knobs[8].state = uninitialized;
  knobs[8].label = "resonance";

  knobs[9].data = 11; /* expression */
  knobs[9].action = &knobVolumeAction;
  knobs[9].state = uninitialized;
  knobs[9].label = "volume";

  knobs[6].data = 91; /* reverb send */
  knobs[6].action = &knobReverbAction;
  knobs[6].state = uninitialized;
  knobs[6].label = "reverb";

  knobs[0].data = 0;
  knobs[0].action = &knobVolumeAttackRateAction;
  knobs[0].state = active;
  knobs[0].label = "volume attack";

  knobs[1].data = 0;
  knobs[1].action = &knobVolumeReleaseRateAction;
  knobs[1].state = active;
  knobs[1].label = "volume release";

  knobs[4].data = 0;
  knobs[4].action = &knobFilterAttackRateAction;
  knobs[4].state = active;
  knobs[4].label = "filter attack";

  knobs[5].data = 0;
  knobs[5].action = &knobFilterReleaseRateAction;
  knobs[5].state = active;
  knobs[5].label = "filter release";

  knobs[2].data = 0;
  knobs[2].action = &knobPressureExponentAction;
  knobs[2].state = active;
  knobs[2].label = "pressure curve";
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
    key->pitch.type = ratioType; \
    key->pitch.ratio.a = a; \
    key->pitch.ratio.b = b; \
    key->originalPitch = key->pitch; \
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
  int pbThresholdPressure = thresholdPressure * 4.0;
  int pbMaxPressure = maxPressure + (pbThresholdPressure - thresholdPressure);   //(thresholdPressure + maxPressure) - (pbThresholdPressure + 500);

  controls[8+16+6][3].update = pbUpUpdate;
  controls[8+16+6][3].thresholdPressure = pbThresholdPressure;
  controls[8+16+6][3].maxPressure = pbMaxPressure;
  controls[8+16+6][3].area = 4.0f;
  controls[8+16+7][3].update = pbDownUpdate;
  controls[8+16+7][3].thresholdPressure = pbThresholdPressure;
  controls[8+16+7][3].maxPressure = pbMaxPressure;
  controls[8+16+7][3].area = 4.0f;
}

void ratioSwap(uint32_t prime, struct RatioPitch to) {
  ratioRestore();
  Serial.println ("replacing prime " + String(prime) + " with " + String(to.a) + "/" + String(to.b));
  for (int i = 0; i < maxKeys; i++) {
    struct Key *key = &keys[i];
    uint32_t a = key->pitch.ratio.a;
    uint32_t b = key->pitch.ratio.b;

    if (key->originalPitch.ratio.a % prime == 0) {
      key->pitch.ratio.a = (a / prime) * to.a;
      key->pitch.ratio.b *= to.b;
    }
    if (key->originalPitch.ratio.b % prime == 0) {
      key->pitch.ratio.b = (b / prime) * to.a;
      key->pitch.ratio.a *= to.b;
    }

    key->pitch.simplify();
  }
}

void ratioRestore() {
  for (int i = 0; i < maxKeys; i++) {
    struct Key *key = &keys[i];
    key->pitch = key->originalPitch;
  }
}

enum substitution {
  subDefault,
  sub7_11,
  sub7_13,
  sub7_25,
  swap5_7
};

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

void versionAction(void *data) {
  status1TextUpdate(fwversion, 2000000);
}

bool enableVisualizer = false;
bool lastEnableVisualizer = enableVisualizer;
bool debugShowResistances = false;
bool debugShowCalibration = false;

float reverbSize = 1.0;
float reverbHiDamp = 0.3;
float reverbLoDamp = 0.1;
float reverbLowPass = 0.2;
float reverbDiffusion = 1.0;

struct MenuItem allNotesOffSlowMenuItem("notes off", allNotesOffSlowAction);
struct MenuItem useUsbMenuItem("usb midi", &useUsbMidi);
struct MenuItem useDinMenuItem("din5 midi", &useDinMidi);
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

struct MenuItem doVelocityMenuItem("velocity", &doMpeDynamicVelocity);
struct MenuItem doPressureMenuItem("pressure", &doMpeDynamicPressure);
struct MenuItem doVelocityMenuItemTerse("vel", &doMpeDynamicVelocity);
struct MenuItem doPressureMenuItemTerse("pre", &doMpeDynamicPressure);
struct MenuItem doPolyAfterTouchMenuItem("poly at", &doMpePolyAfterTouch);
struct MenuItem pressureBackoffMenuItem("p backoff", value, &pressureBackoff);
struct MenuItem bendUpOnlyMenuItem("upbend only", &bendUpOnly);
struct MenuItem bendDownOnlyMenuItem("dnbend only", &bendDownOnly);
struct MenuItem enableBenderMenuItem("bender", &enableBender);
struct MenuItem lockMenuItem("lock", &lock);

struct MenuItem mpeBankLsbMenuItem("bank LSB", value, &mpeBankLsb, &mpeBankLsbMin, &mpeBankLsbMax);
struct MenuItem mpeBankMsbMenuItem("bank MSB", value, &mpeBankMsb, &mpeBankMsbMin, &mpeBankMsbMax);

struct MenuItem debugShowResistancesMenuItem("show res", &debugShowResistances);
struct MenuItem debugShowCalibrationMenuItem("show cal", &debugShowCalibration);

struct MenuItem unlockBankRangeMenuItem("unlock", &unlockBankRange);
struct MenuItem sendETuningTableMenuItem("tx E! tt", sendETuningTableAction);

struct MenuItem doCCPassThroughMenuItem("CC thru", &doCCPassThrough);

struct MenuItem versionMenuItem("fw version", versionAction);

uint32_t zero = 0;
uint32_t one = 1;
uint32_t hundred = 100;
uint32_t maxMidiValue = 127;
uint32_t maxChannels = 16;
uint32_t substitutions = subDefault;
uint32_t substitutionsActive = subDefault;

#include <synth_waveform.h>

//int synthWaveform = WAVEFORM_SAWTOOTH;
int synthWaveform = WAVEFORM_TRIANGLE_VARIABLE;


struct MenuItem pressureCCMenuItem("pressure CC", value, &pressureCC, &zero, &maxMidiValue);
struct MenuItem mpeProgramChangeMenuItem("patch", value, &programChange, &zero, &maxMidiValue);
struct MenuItem mpeChannelsMenuItem("channels", value, &mpeChannels, &zero, &maxChannels);
struct MenuItem minVelocityMenuItem("min vel", value, &minVelocity, &one, &maxMidiValue);

struct MenuItem surgeXtPresetMenuItem("Surge XT", selection, &mpeSettings, (uint32_t)&mpeSettingsSurgeXT);
struct MenuItem kspPresetMenuItem("Keystep Pro", selection, &mpeSettings, (uint32_t)&mpeSettingsKSP);
struct MenuItem rd300PresetMenuItem("RD-300NX", selection, &mpeSettings, (uint32_t)&mpeSettingsRD300NX);
struct MenuItem xv2020PresetMenuItem("XV-2020", selection, &mpeSettings, (uint32_t)&mpeSettingsXV2020);
struct MenuItem sp300PresetMenuItem("SP-300", selection, &mpeSettings, (uint32_t)&mpeSettingsSP300);
struct MenuItem trinityPresetMenuItem("Trinity", selection, &mpeSettings, (uint32_t)&mpeSettingsTrinity);
struct MenuItem x50PresetMenuItem("X50", selection, &mpeSettings, (uint32_t)&mpeSettingsX50);
struct MenuItem fb01PresetMenuItem("FB-01", selection, &mpeSettings, (uint32_t)&mpeSettingsFB01);
struct MenuItem dx7EPresetMenuItem("DX7 with E!", selection, &mpeSettings, (uint32_t)&mpeSettingsDx7E);
struct MenuItem moxPresetMenuItem("MOX6/8", selection, &mpeSettings, (uint32_t)&mpeSettingsMox8);
struct MenuItem proteus2000PresetMenuItem("Proteus 2k", selection, &mpeSettings, (uint32_t)&mpeSettingsProteus);
struct MenuItem phattyPresetMenuItem("Slim Phatty", selection, &mpeSettings, (uint32_t)&mpeSettingsPhatty);
struct MenuItem dexedPresetMenuItem("Dexed", selection, &mpeSettings, (uint32_t)&mpeSettingsDexed);
struct MenuItem qSynthPresetMenuItem("QSynth", selection, &mpeSettings, (uint32_t)&mpeSettingsQSynth);
struct MenuItem mt2PresetMenuItem("MT v2", selection, &mpeSettings, (uint32_t)&mpeSettingsMT2);

struct MenuItem arturiaMenu("Arturia", submenu, &kspPresetMenuItem);
struct MenuItem befacoMenu("Befaco", submenu, &mt2PresetMenuItem);
struct MenuItem emuMenu("E-mu", submenu, &proteus2000PresetMenuItem);
struct MenuItem korgMenu("Korg", submenu, &sp300PresetMenuItem, &trinityPresetMenuItem, &x50PresetMenuItem);
struct MenuItem moogMenu("Moog", submenu, &phattyPresetMenuItem);
struct MenuItem rolandMenu("Roland", submenu, &rd300PresetMenuItem, &xv2020PresetMenuItem);
struct MenuItem yamahaMenu("Yamaha", submenu, &dx7EPresetMenuItem, &fb01PresetMenuItem, &moxPresetMenuItem);

struct MenuItem* brandsMenu[] = {&arturiaMenu, &befacoMenu, &dexedPresetMenuItem, &emuMenu, &korgMenu, &moogMenu, &rolandMenu, &surgeXtPresetMenuItem, &qSynthPresetMenuItem, &yamahaMenu};
struct MenuItem outputPresetsMenu("dev presets", submenu, &brandsMenu[0], 10);
struct MenuItem pbRangeMenu("bend range", submenu, &pb2MenuItem, &pb7MenuItem, &pb12MenuItem, &pb24MenuItem, &pb48MenuItem);
struct MenuItem noteOnFirstMenuItem("note-on 1st", toggle, &noteOnFirst);
struct MenuItem maxPressureMenuItem("max p", value, &maxMpePressure, &zero, &maxMidiValue);
struct MenuItem minPressureMenuItem("min p", value, &minMpePressure, &zero, &maxMidiValue);

struct MenuItem* outputMenuItems[] = {&useUsbMenuItem, &useDinMenuItem, &outputPresetsMenu, &mpeHandshakeMenuItem, &mpeChannelsMenuItem, &pbRangeMenu, &pressureCCMenuItem, &maxPressureMenuItem, &minPressureMenuItem, &minVelocityMenuItem, &doCCPassThroughMenuItem};
struct MenuItem outputMenu("midi", submenu, &outputMenuItems[0], 11);

struct MenuItem sub7_11MenuItem("7/4->11/8", selection, &substitutions, sub7_11);
struct MenuItem swap5_7MenuItem("5:7 swap", selection, &substitutions, swap5_7);
struct MenuItem sub7_13MenuItem("14/9->13/8", selection, &substitutions, sub7_13);
struct MenuItem sub7_25MenuItem("7/4->25/16", selection, &substitutions, sub7_25);
struct MenuItem subDefaultMenuItem("default", selection, &substitutions, subDefault);
struct MenuItem ratioSubsMenu("substitute", submenu, &sub7_11MenuItem, &sub7_13MenuItem, &sub7_25MenuItem, &swap5_7MenuItem, &subDefaultMenuItem);

struct MenuItem* controlsMenuItems[] = {&doVelocityMenuItem, &doPressureMenuItem, &doPolyAfterTouchMenuItem, &pressureBackoffMenuItem, &enableBenderMenuItem, &ratioSubsMenu};
struct MenuItem controlsMenu("controls", submenu, &controlsMenuItems[0], 6);

struct MenuItem visualizerMenuItem("visualizer", toggle, &enableVisualizer);
struct MenuItem screenBrightnessMenu("brightness", submenu, &screen10MenuItem, &screen25MenuItem, &screen50MenuItem, &screen75MenuItem, &screen100MenuItem);
struct MenuItem interfaceMenu("interface", submenu, &screenBrightnessMenu, &visualizerMenuItem);
struct MenuItem patchesMenu("patches", submenu, &mpeBankMsbMenuItem,  &mpeBankLsbMenuItem, &mpeProgramChangeMenuItem, &unlockBankRangeMenuItem);

struct MenuItem sawMenuItem("sawtooth", selection, &synthWaveform, WAVEFORM_BANDLIMIT_SAWTOOTH);
struct MenuItem trileanMenuItem("skew", selection, &synthWaveform, WAVEFORM_TRIANGLE_VARIABLE);
struct MenuItem triMenuItem("triangle", selection, &synthWaveform, WAVEFORM_TRIANGLE);
struct MenuItem pulseMenuItem("pulse", selection, &synthWaveform, WAVEFORM_BANDLIMIT_PULSE);
struct MenuItem squareMenuItem("square", selection, &synthWaveform, WAVEFORM_BANDLIMIT_SQUARE);
struct MenuItem sineMenuItem("sine", selection, &synthWaveform, WAVEFORM_SINE);

struct MenuItem reverbSizeMenuItem("size", &reverbSize);
struct MenuItem reverbHiDampMenuItem("high damp", &reverbHiDamp);
struct MenuItem reverbLoDampMenuItem("low damp", &reverbLoDamp);
struct MenuItem reverbLowPassMenuItem("low pass", &reverbLowPass);
struct MenuItem reverbDiffusionMenuItem("diffusion", &reverbDiffusion);

struct MenuItem* waveformMenuItems[] = {&sawMenuItem, &trileanMenuItem, &triMenuItem, &sineMenuItem, &squareMenuItem, &pulseMenuItem};

struct MenuItem oscillatorMenu("oscillator", submenu, waveformMenuItems, 6);
struct MenuItem reverbMenu("reverb", submenu, &reverbSizeMenuItem, &reverbHiDampMenuItem, &reverbLoDampMenuItem, &reverbLowPassMenuItem, &reverbDiffusionMenuItem);

struct MenuItem synthMenu("synth", submenu, &oscillatorMenu, &reverbMenu);

struct MenuItem* debugMenuItems[] = {&versionMenuItem, &debugShowResistancesMenuItem, &debugShowCalibrationMenuItem, &noteOnFirstMenuItem, &bendUpOnlyMenuItem, &bendDownOnlyMenuItem, &lockMenuItem};
struct MenuItem debugMenu("debug", submenu, &debugMenuItems[0], 7);

struct MenuItem configMenu("settings", submenu, &outputMenu, &controlsMenu, &interfaceMenu, &synthMenu, &debugMenu);

struct MenuItem rootMenu("", submenu, &configMenu, &patchesMenu, &emptyMenuItem, &emptyMenuItem, &allNotesOffSlowMenuItem);

void menuSetup() {
  menuSelect(&rootMenu, 0);

  patchesMenu.addOption(&doVelocityMenuItemTerse, 0);
  patchesMenu.addOption(&doPressureMenuItemTerse, 2);
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

  /* no shift register output bit set yet */
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
          redrawVisualizer();
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

  synthUpdate(delta);
  mpeUpdate(delta);

  if (substitutions != substitutionsActive) {
    struct RatioPitch ratio;
    switch (substitutions) {
      case subDefault:
        ratioRestore();
        break;
      case sub7_11: /* convert 7/4 into 11/4 */
        ratio.a = 11;
        ratio.b = 1;
        ratioSwap(7, ratio);
        break;
      case sub7_13: /* convert 14,9 into 13,8 */
        ratio.a = 13*9;
        ratio.b = 16;
        ratioSwap(7, ratio);
        break;
      case sub7_25:
        ratio.a = 25*4;
        ratio.b = 16;
        ratioSwap(7, ratio);
        break;
      case swap5_7:
        ratio.a = 36;
        ratio.b = 5;
        ratioSwap(7, ratio);
        ratio.a = 36;
        ratio.b = 7;
        ratioSwap(5, ratio);
      default:
        Serial.println("unexpected substitution type");
        break;
    }
    substitutionsActive = substitutions;
  }

  if (verbose  && debugShowResistances) {
    showResistances();
  }

  if (iteration % 16 == 0) { 
    int surplus = noteOnCount - noteOffCount;
    if (surplus > 10) {
      surplus = 5;
    }

    /* todo: only set if the value changes */
    setLed(5, surplus);
    showLeds();
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
        if (doCCPassThrough) {
          mpeMulticastCC(cc, data);
        }
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

  /*
  if (verbose && true) {
    Serial.println("MIDI channel state:");
    for (int i=0; i < 16; i++) {
      struct MpeChannelState *chan = &mpeState[i];
      Serial.println("playing: " + String(chan->playing) + " vol: " + String(chan->volume) + " age: " + String(chan->age));
    }
  } */

  if (verbose) {
    Serial.println(String("MIDI Sent: ") + midiMsgsSent + " received: " + midiMsgsReceived + " vol release rate " + String(volumeReleaseRate) + " pressure exponent " + String(pressureExponent) + " pb +/- " + String(pbUp) + "  " + String(pbDown));
  }
  iteration++;
}
