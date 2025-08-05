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

void setLed(uint16_t led, uint32_t color) {
  leds.setPixel(led, color);
}

void showLeds() {
  leds.show();
}
