/*
Copyright 2025 Jim Snow, Desiderata Systems LLC

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

/* AUDIO */

#include <Audio.h>
#include <Wire.h>
#include <SPI.h>
#include <SD.h>
#include <SerialFlash.h>
#include "tapered_synth_waveform.h"
#include "effect_platervbstereo.h"

#define synthVoices 16
#define mixers ((synthVoices+3)/4)  /* each of the first-level mixers take up to 4 inputs */

// GUItool: begin automatically generated code
TaperedAudioSynthWaveform       waveform[synthVoices];      //xy=106.33336639404297,53.3333683013916
AudioSynthKarplusStrong  string[synthVoices];
AudioFilterBiquad        biquad[synthVoices];        //xy=263.33333587646484,54.33335304260254
AudioMixer4              mixer[mixers];         //xy=418.33334159851074,74.33335876464844
AudioMixer4              perVoiceMixer[synthVoices];
AudioSynthWaveformSine   sine1;          //xy=419.3333435058594,285.33337020874023
AudioMixer4              mixer4;         //xy=611.333381652832,178.33335494995117
//AudioEffectReverb        reverb1;        //xy=707.3333854675293,61.3333854675293
//AudioEffectFreeverb      reverb1;
AudioEffectPlateReverb   reverb1;
AudioMixer4              mixerL;         //xy=841.3333892822266,185.33335494995117
AudioMixer4              mixerR;
AudioOutputI2S2          i2s2_1;         //xy=1014.3333892822266,243.33335876464844

AudioConnection          patchCord1[synthVoices] = {};
AudioConnection          patchCord2[synthVoices] = {};
AudioConnection          patchCord3[mixers];
AudioConnection          patchCord6(sine1, 0, mixerL, 2);
AudioConnection          patchCord7(sine1, 0, mixerR, 2);

AudioConnection          perVoiceMixerPatchCord[synthVoices] = {};
AudioConnection          stringPatchCord[synthVoices] = {};

//AudioConnection          pc1(waveform[0], 0, mixer[3], 0);
//AudioConnection          pc2(waveform[1], 0, mixer[3], 1);
//AudioConnection          pc3(waveform[2], 0, mixer[3], 2);
//AudioConnection          pc4(waveform[3], 0, mixer[3], 3);

//AudioConnection          pcmm1(mixer[0], 0, mixer4, 0);
//AudioConnection          pcmm2(mixer[1], 0, mixer4, 1);
//AudioConnection          pcmm3(mixer[2], 0, mixer4, 2);
//AudioConnection          pcmm4(mixer[3], 0, mixer4, 3);

//AudioConnection          patchCord1(waveform1, biquad1);
//AudioConnection          patchCord2(biquad1, 0, mixer1, 0);
//AudioConnection          patchCord3(mixer1, 0, mixer4, 0);
//AudioConnection          patchCord4(mixer3, 0, mixer4, 2);
//AudioConnection          patchCord5(mixer2, 0, mixer4, 1);
//AudioConnection          patchCord6(sine1, 0, mixer4, 3);

AudioConnection          patchCord8(mixer4, reverb1);
AudioConnection          patchCord9(mixer4, 0, mixerL, 1);
AudioConnection          patchCord10(mixer4, 0, mixerR, 1);
AudioConnection          patchCord11(reverb1, 0, mixerL, 0);
AudioConnection          patchCord12(reverb1, 1, mixerR, 0);
AudioConnection          patchCord13(mixerL, 0, i2s2_1, 0);
AudioConnection          patchCord14(mixerR, 0, i2s2_1, 1);
AudioControlSGTL5000     sgtl5000_1;     //xy=969.3333435058594,44.33335590362549
// GUItool: end automatically generated code

/* ON-BOARD SYNTHESIZER */


float keyFilterWeight = 0.2f;
float filterPitchTracking = 0.5f;

struct SynthVoice {
  uint16_t owner;
  TaperedAudioSynthWaveform *osc;
  AudioFilterBiquad *filter;
  AudioMixer4 *mixer;
  int mixerChannel;
  uint32_t age;
  float volume;
  float pitch; // relative to 1:1
  float frequency;  // in hz, as of when the note began
  float filterAmount;
  AudioSynthKarplusStrong *string;
};

struct SynthVoice voices[synthVoices] = {};

float synthOscillatorMaxVolume = 0.2f;
float masterSynthVolume = 1.0f;

float audioTaper(float amount) {
  amount = clamp(amount);
  return pow(amount, 3.0f);
}

void setSynthVolume(float amount) {
  masterSynthVolume = reverseLerp(0.0f, amount*amount, 0.5f);

  AudioNoInterrupts();
  for (int p = 0; p < 4; p++) {
    mixer4.gain(p, masterSynthVolume * 4.0f);
  }
  AudioInterrupts();
}

void setSynthReverb(float amount) {
  if (amount > 1.0f) {
    amount = 1.0f;
  }

  if (amount < 0.0f) {
    amount = 0.0f;
  }

  amount = audioTaper(amount);

  float dry = 1.0f - amount;
  float wet = amount * 3.0;  /* plate reverb is a little quiet */

  AudioNoInterrupts();
  mixerL.gain(0, wet);
  mixerR.gain(0, wet);
  mixerL.gain(1, dry);
  mixerR.gain(1, dry);
  AudioInterrupts();
}

void setVoiceFilter(struct SynthVoice* voice, float keyAmount, float knobAmount, float res) {
  knobAmount = audioTaper(knobAmount);

  float amount = ((knobAmount * (1.0 + keyFilterWeight) - keyFilterWeight) * (1.0f - keyFilterWeight)) + (keyAmount * keyFilterWeight);
  amount *= (1.0f - filterPitchTracking) + (voice->pitch * filterPitchTracking);

  float freq = reverseLerp(20.0f, amount, 20000.0f);
  res = reverseLerp(0.1f, res, 5.0f);
  AudioNoInterrupts();
  voice->filter->setLowpass(0, freq, res);
  AudioInterrupts();
  //Serial.println("set filter freq " + String(freq) + " res " + String(res));
}

void setSynthPitchBend(float pitch) {
  AudioNoInterrupts();
  for (int v = 0; v < synthVoices; v++) {
    struct SynthVoice *voice = &voices[v];
    if (voice->volume > 0.0f) {
      voice->osc->frequency(pitch * voice->frequency);
    }
  }
  AudioInterrupts();
}

float modRate = 240.0; /* bpm */
float modDepth = 0.0f;
float modValue = 1.0f;

void setSynthModulation(float mod) {
  modDepth = clamp(mod);
}

/* m=minutes elapsed, return a value between 0 and 1 */
float modFunction(float m) {
  return (sin(m * modRate * 2.0f * M_PI) + 1.0f) * 0.5f;
}

void updateModulationValue() {
  float result = modFunction((float)micros() / (1000000.0f * 60.0f));
  modValue = 1.0f - (result * modDepth);
}

void setVoiceVolume(struct SynthVoice *voice, float volume, uint32_t deltaUsecs) {
  volume = clamp(volume);
  volume = pow(volume, pressureExponent);

  float maxVolume = voice->volume + volumeAttackRate  * ((float)deltaUsecs/1000000.0f) * (volume * volume);
  if (volume > maxVolume) {
    volume = maxVolume;
  }

  float minVolume = voice->volume - (volumeReleaseRate * (voice->volume + 0.2f)) * ((float)deltaUsecs/1000000.0f);
  if (minVolume < 0.0f) {
    minVolume = 0.0f;
  }

  if (volume < minVolume) {
    volume = minVolume;
  }

  AudioNoInterrupts();
  voice->osc->amplitude(audioTaper(volume) * synthOscillatorMaxVolume * modValue);
  voice->volume = volume;
  AudioInterrupts();
}

void setVoiceFilterMod(struct SynthVoice *voice, float pressure, uint32_t deltaUsecs) {
  pressure = clamp(pressure);
  pressure = pow(pressure, pressureExponent);

  float minFilter = voice->filterAmount - (filterReleaseRate * (voice->filterAmount + 0.2f)) * ((float)deltaUsecs/1000000.0f);
  float maxFilter = voice->filterAmount + filterAttackRate  * ((float)deltaUsecs/1000000.0f) * (pressure);
  float filter;

  if (pressure < minFilter) {
    filter = minFilter;
  } else if (pressure > maxFilter) {
    filter = maxFilter;
  } else {
    filter = pressure;
  }
  voice->filterAmount = filter;
  setVoiceFilter(voice, filter, filterAmount, resonanceAmount);
}

struct SynthVoice *getSynthVoice(uint16_t owner) {
  int best = 0;
  float bestVolume = 1.0f;
  for (int v = 0; v < synthVoices; v++) {
    Serial.println("getSynthVoice voice " + String(v) + " volume " + String(voices[v].volume));
    if (voices[v].volume < bestVolume) {
      best = v;
      bestVolume = voices[v].volume;
      Serial.println("bestVolume is " + String(bestVolume));
    }
  }
  voices[best].owner = owner;

  Serial.println("most idle synth voice " + String(best) + " volume " + String(bestVolume));
  return &voices[best];
}

struct SynthVoice *beginSynthNote(double pitch, float velocity, float pressure, uint16_t owner) {
  struct SynthVoice* voice = getSynthVoice(owner);

  float bend = calculatePitchBend(pbUp, pbDown);
  voice->pitch = pitch;
  voice->frequency = pitch * pitchReferenceHz();
  voice->osc->frequency(voice->frequency * bend);

  AudioNoInterrupts();
  voice->osc->begin(synthWaveform);
  AudioInterrupts();

  setVoiceVolume(voice, pressure, 100);
  setVoiceFilterMod(voice, pressure, 100);
  //voice->string->noteOn(voice->pitch * bend * pitchReferenceHz(), audioTaper(velocity));

  return voice;
}

void continueSynthNote(struct SynthVoice* voice, float pressure, uint32_t deltaUsecs, uint16_t owner) {
  if (voice->owner != owner) {
    return;
  }
  //Serial.println("setting synth volume " + String(pressure));
  setVoiceVolume(voice, pressure, deltaUsecs);
  setVoiceFilterMod(voice, pressure, deltaUsecs);
  return;
}

void endSynthNote(struct SynthVoice* voice, uint16_t owner) {
  if (voice->owner != owner) {
    return;
  }
  voice->owner = noOne;
  //voice->string->noteOff(0.5f);
}

void synthSetup() {
  //patchCord1[0] = AudioConnection(waveform[0], 0, mixer5, 3);

  for (int v=0; v<synthVoices; v++) {
    int m = v/4;
    int port = v%4;

    patchCord1[v].connect(waveform[v], biquad[v]);
    patchCord2[v].connect(biquad[v], 0, perVoiceMixer[v], 0);
    perVoiceMixer[v].gain(0, 1.0f);
    stringPatchCord[v].connect(string[v], 0, perVoiceMixer[v], 1);
    perVoiceMixer[v].gain(1, 0.3f);
    perVoiceMixerPatchCord[v].connect(perVoiceMixer[v], 0, mixer[m], port);

    //patchCord2[v].connect(biquad[v], 0, mixer[m], port);

    //Serial.println("connecting waveform " + String(v) + " to mixer " + String(m) + " port " + String(port));

    voices[v].osc = &waveform[v];
    voices[v].filter = &biquad[v];
    voices[v].mixer = &mixer[m];
    voices[v].mixerChannel = v%4;
    voices[v].age = 0xffffffff;
    voices[v].volume = 0.0f;
    voices[v].string = &string[v];

    voices[v].osc->begin(0.0f, 440.0f, WAVEFORM_SAWTOOTH);
    voices[v].filter->setLowpass(0, 800.0f, 0.1f);
  }


  // doesn't actually affect volume
  sgtl5000_1.volume(1.0);
    
  //reverb1.roomsize(1.0f);
  //reverb1.damping(0.9f);

  reverb1.size(reverbSize);
  reverb1.lowpass(reverbLowPass);
  reverb1.lodamp(reverbLoDamp);
  reverb1.hidamp(reverbHiDamp);
  reverb1.diffusion(reverbDiffusion);

  for (int m=0; m < mixers; m++) {
    patchCord3[m].connect(mixer[m], 0, mixer4, m);
  }

  setSynthVolume(0.5f);
  setSynthReverb(0.4f);
}

float lastReverbSize = reverbSize;
float lastReverbHiDamp = reverbHiDamp;
float lastReverbLoDamp = reverbLoDamp;
float lastReverbLowPass = reverbLowPass;
float lastReverbDiffusion = reverbDiffusion;

void synthUpdate(uint32_t deltaUsecs) {
  if (!doLocalSynth) {
    return;
  }

  setSynthPitchBend(calculatePitchBend(pbUp, pbDown));
  updateModulationValue();
  AudioNoInterrupts();
  for (int v=0; v < synthVoices; v++) {
    struct SynthVoice *voice = &voices[v];

    if (voice->owner == noOne) {
      setVoiceVolume(voice, 0.0f, deltaUsecs);
    }
  }

  if (reverbSize != lastReverbSize) {
    reverb1.size(reverbSize);
    lastReverbSize = reverbSize;
  }

  if (reverbHiDamp != lastReverbHiDamp) {
    reverb1.hidamp(reverbHiDamp);
    lastReverbHiDamp = reverbHiDamp;
  }

  if (reverbLoDamp != lastReverbLoDamp) {
    reverb1.lodamp(reverbLoDamp);
    lastReverbLoDamp = reverbLoDamp;
  }

  if (reverbLowPass != lastReverbLowPass) {
    reverb1.lowpass(reverbLowPass);
    lastReverbLowPass = reverbLowPass;
  }

  if (reverbDiffusion != lastReverbDiffusion) {
    reverb1.diffusion(reverbDiffusion);
    lastReverbDiffusion = reverbDiffusion;
  }
  AudioInterrupts();
}

void audioSetup() {
  AudioMemory(40);
  AudioNoInterrupts();
  synthSetup();
  /* a440 test tone */
  sine1.frequency(440.0f);
  sine1.amplitude(0.0f);
  AudioInterrupts();

  //pinMode(mutePin, OUTPUT);
  //digitalWrite(mutePin, HIGH); /* unmute */
  Serial.println("audio setup complete");
}

void setPitchReference(double freq) {
  AudioNoInterrupts();
  sine1.frequency(freq);
  sine1.amplitude(0.0f);
  AudioInterrupts();
}

