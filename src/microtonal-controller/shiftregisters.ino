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

void shiftRegisterSetup() {
  pinMode(shiftRegisterOutPin, OUTPUT);
  pinMode(shiftRegisterClockPin, OUTPUT);
  shiftRegisterReset(0);
}

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
