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

/* ADCs */

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


