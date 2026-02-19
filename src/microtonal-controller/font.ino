#include <ILI9341_t3.h>

const int justifyLeft = 0;
const int justifyCenter = 1;
const int justifyRight = 2;

struct CharGlyph {
  CharGlyph(uint8_t width,
            uint8_t fg0, uint8_t mid0, 
            uint8_t fg1, uint8_t mid1,
            uint8_t fg2, uint8_t mid2,
            uint8_t fg3, uint8_t mid3,
            uint8_t fg4, uint8_t mid4,
            uint8_t fg5, uint8_t mid5,
            uint8_t fg6, uint8_t mid6,
            uint8_t fg7, uint8_t mid7,
            uint8_t fg8, uint8_t mid8,
            uint8_t fg9 = 0, uint8_t mid9 = 0,
            uint8_t fg10 = 0, uint8_t mid10 = 0,
            uint8_t fg11 = 0, uint8_t mid11 = 0,
            uint8_t fg12 = 0, uint8_t mid12 = 0) : width{width}, shift{width} {
    data[0] = fg0;
    data[1] = mid0;
    data[2] = fg1;
    data[3] = mid1;
    data[4] = fg2;
    data[5] = mid2;
    data[6] = fg3;
    data[7] = mid3;
    data[8] = fg4;
    data[9] = mid4;
    data[10] = fg5;
    data[11] = mid5;
    data[12] = fg6;
    data[13] = mid6;
    data[14] = fg7;
    data[15] = mid7;
    data[16] = fg8;
    data[17] = mid8;
    data[18] = fg9;
    data[19] = mid9;
    data[20] = fg10;
    data[21] = mid10;
    data[22] = fg11;
    data[23] = mid11;
    data[24] = fg12;
    data[25] = mid12;
  }

  CharGlyph(uint8_t width = 0) : width{width}, shift{width} {
    for (int i=0; i < 26; i++) {
      data[i] = 0;
    }
  }

  bool isFg(uint16_t x, uint16_t y) {
    if (x >= width) {
      if (next == nullptr) {
        return false;
      } else {
        return next->isFg(x-shift, y);
      }
    }

    uint8_t row = data[y * 2];
    int bitIndex = width - (x + 1);
    if (bitIndex < 0 || bitIndex > 7) {
      return false;
    }
    return (row & (1 << bitIndex)) > 0;
  }

  bool isMid(uint16_t x, uint16_t y) {
    if (x >= width) {
      if (next == nullptr) {
        return false;
      } else {
        return next->isMid(x-shift, y);
      }
    }
    uint8_t row = data[y * 2 + 1];
    int bitIndex = width - (x + 1);
    if (bitIndex < 0 || bitIndex > 7) {
      return false;
    }
    return (row & (1 << bitIndex)) > 0;
  }

  /* only flips top 9 pixels */
  void flipUpDown() {
    for (int i=0; i<4; i++) {
      uint8_t a1 = data[i * 2];
      uint8_t a2 = data[i * 2 + 1];
      
      data[i * 2] = data[(8 - i) * 2];
      data[i * 2 + 1] = data[(8 - i) * 2 + 1];

      data[(8 - i) * 2] = a1;
      data[(8 - i) * 2 + 1] = a2;
    }
  }
  void flipRL() {

  }

  int fullWidth() {
    CharGlyph *g = this;
    int fw = 0;
    do {
      fw += g->width;
      g = g->next;
    } while (g != nullptr);

    return fw;
  }

  int fullShift() {
    CharGlyph *g = this;
    int fs = 0;
    do {
      fs += g->shift;
      g = g->next;
    } while (g != nullptr);

    return fs;
  }

  uint8_t data[26];
  uint8_t width;  /* how many pixels wide data is */
  uint8_t shift;  /* how many pixels to shift over when printing */
  CharGlyph* next = nullptr; /* chain glyphs if you need more than 8 pixels */
};

struct CharGlyph glyphs[128];
struct CharGlyph wExtra;
struct CharGlyph atExtra;

struct CharGlyph din5PlugSmall(8,
  0b00111110, 0b01000001,
  0b01011101, 0b10000000,
  0b10000000, 0b01000001,
  0b10000000, 0b00000000,
  0b10100010, 0b01000001,
  0b10000000, 0b00010100,
  0b10101010, 0b01000001,
  0b01000001, 0b10101010,
  0b00111110, 0b01000001);

struct CharGlyph din5PlugSmall2(1,
  0b0, 0b0,
  0b0, 0b1,
  0b1, 0b0,
  0b1, 0b0,
  0b1, 0b0,
  0b1, 0b0,
  0b1, 0b0,
  0b0, 0b1,
  0b0, 0b0);

struct CharGlyph din5Plug(8,
  0b00001110, 0b00010001,
  0b00111111, 0b01000000,
  0b01001110, 0b00100000,
  0b01000000, 0b10000000,
  0b10000000, 0b00010001,
  0b10100000, 0b00000000,
  0b10000000, 0b01001010,
  0b01010001, 0b10000100,
  0b01000100, 0b00100000,
  0b00110001, 0b01000100,
  0b00001110, 0b00010001);

struct CharGlyph din5Plug2(3,
  0b000, 0b000,
  0b100, 0b010,
  0b010, 0b100,
  0b010, 0b001,
  0b001, 0b000,
  0b101, 0b000,
  0b001, 0b010,
  0b010, 0b001,
  0b010, 0b100,
  0b100, 0b010,
  0b000, 0b000);

struct CharGlyph usbSymbol(6,
  0b000000, 0b0,
  0b000000, 0b0,
  0b011100, 0b00000000,
  0b111110, 0b10001000,
  0b111111, 0b00000000,
  0b111110, 0b00000000,
  0b011100, 0b10001000,
  0b000000, 0b0,
  0b000000, 0b0);

struct CharGlyph usbSymbol2(8,
  0b00000000, 0b0,
  0b00011111, 0b00100000,
  0b00100000, 0b01010000,
  0b01000000, 0b10100000,
  0b11111111, 0b00000000,
  0b00001000, 0b00010100,
  0b00000100, 0b00001010,
  0b00000011, 0b00000100,
  0b00000000, 0b0);

struct CharGlyph usbSymbol3(8,
  0b01100000, 0b10010000,
  0b11110000, 0b00000000,
  0b01100000, 0b10010000,
  0b00000000, 0b00000000,
  0b11111111, 0b00000000,
  0b00000000, 0b00000000,
  0b00111100, 0b00000000,
  0b11111100, 0b0,
  0b00111100, 0b0);

struct CharGlyph usbSymbol4(4,
  0b0000, 0b0,
  0b0000, 0b10000000,
  0b1000, 0b01000000,
  0b1100, 0b00100000,
  0b1110, 0b00010000,
  0b1100, 0b00100000,
  0b1000, 0b01000000,
  0b0000, 0b10000000,
  0b0000, 0b0);

struct CharGlyph downSymbol(8,
  0b00000000, 0b00000000,
  0b00001111, 0b00010000,
  0b00011110, 0b00000000,
  0b00011110, 0b00100000,
  0b00110000, 0b00001000,
  0b00111000, 0b01000100,
  0b01111100, 0b00000010,
  0b01111110, 0b10000001,
  0b11111111, 0b00000000);

struct CharGlyph downSymbol2(8,
  0b00000000, 0b00000000,
  0b11110000, 0b00001000,
  0b01111000, 0b00000000,
  0b01111000, 0b00000100,
  0b00001100, 0b00010000,
  0b00011100, 0b00100010,
  0b00111110, 0b01000000,
  0b01111110, 0b10000001,
  0b11111111, 0b00000000);

struct CharGlyph upSymbol = downSymbol;
struct CharGlyph upSymbol2 = downSymbol2;

struct CharGlyph leftSymbol(8,
  0b11111110, 0b00000000,
  0b11111110, 0b00000001,
  0b11110111, 0b00001000,
  0b11100111, 0b00010000,
  0b11000111, 0b00100000,
  0b10000001, 0b01000000,
  0b10000001, 0b01000000,
  0b11000111, 0b00100000,
  0b11100111, 0b00010000,
  0b11110111, 0b00001000,
  0b11111110, 0b00000001,
  0b11111110, 0b00000000);

struct CharGlyph leftSymbol2(2,
  0b00, 0b00,
  0b00, 0b00,
  0b00, 0b00,
  0b00, 0b10,
  0b10, 0b00,
  0b10, 0b01,
  0b10, 0b01,
  0b10, 0b00,
  0b00, 0b10,
  0b00, 0b00,
  0b00, 0b00,
  0b00, 0b00);

struct CharGlyph lockSymbol(8,
  0b00000000, 0b0,
  0b00111100, 0b01000010,
  0b01000010, 0b00100100,
  0b01000010, 0b0,
  0b01000010, 0b0,
  0b11111111, 0b0,
  0b11111111, 0b0,
  0b11111111, 0b0,
  0b11111111, 0b0);

struct CharGlyph leftArrow(8,
  0b0, 0b0,
  0b0, 0b0,
  0b0, 0b0,
  0b00000000, 0b00010000,
  0b00010000, 0b00100000,
  0b00110000, 0b01000000,
  0b01111110, 0b10000001,
  0b00110000, 0b01000000,
  0b00010000, 0b00100000,
  0b00000000, 0b00010000);

struct CharGlyph rightArrow(8,
  0b0, 0b0,
  0b0, 0b0,
  0b0, 0b0,
  0b00000000, 0b00001000,
  0b00001000, 0b00000100,
  0b00001100, 0b00000010,
  0b01111110, 0b10000001,
  0b00001100, 0b00000010,
  0b00001000, 0b00000100,
  0b00000000, 0b00001000);

struct CharGlyph flatSymbol(6,
  0b000000, 0b0,
  0b100000, 0b0,
  0b100000, 0b0,
  0b100000, 0b0,
  0b100000, 0b0,
  0b100000, 0b000000,
  0b101110, 0b010001,
  0b110001, 0b001010,
  0b100010, 0b000101,
  0b111100, 0b000010);

struct CharGlyph centSymbol(7,
  0b0, 0b0,
  0b0, 0b0,
  0b0001000, 0b0,
  0b0111110, 0b1000001,
  0b1001001, 0b0100010,
  0b1001000, 0b0000000,
  0b1001000, 0b0000000,
  0b1001001, 0b0100010,
  0b0111110, 0b1000001,
  0b0001000, 0b0000000);


int charSpacing = 2;
int fontLookBack = 5; /* If a font has pixels up to this far to the left of where "width" implies, we draw those too. */

int drawGlyph(ILI9341_t3 *tft, struct CharGlyph *g, int xOrigin, int yOrigin, uint16_t fgColor, uint16_t midColor, uint16_t bgColor, bool flipV = false, bool flipH = false) {
  int offset = 0;
  int width = g->fullWidth();
  do {
    uint8_t shift = g->fullShift();
    for (int y = 0; y < 12; y++) {
      for (int x = -fontLookBack; x < width; x++) {
        uint16_t color = bgColor;
        int x_ = flipH ? width - (x + 1) : x;
        int y_ = flipV ? 11 - y: y;
        if (g->isFg(x_, y_)) {
          color = fgColor;
        } else if (g->isMid(x_, y_)) {
          color = midColor;
        } else {
          continue;
        }

        tft->drawPixel(x + xOrigin, y + yOrigin, color);
      }
    }

    //g = g->next;
    xOrigin += shift;
    offset += shift;

  } while (false && g != nullptr);

  return offset;
}

struct CharGlyph* glyphLookup(const char *s, int *i, bool *vFlip = nullptr, bool* hFlip = nullptr) {
  CharGlyph* g;
  if (s[*i] == '\\') {
    if (strncmp(&s[*i+1], "usb", 3) == 0) {
      g = &usbSymbol;
      *i += 3;
    } else if(strncmp(&s[*i+1], "din5", 4) == 0) {
      g = &din5Plug;
      *i += 4;
    } else if(strncmp(&s[*i+1], "upkey", 5) == 0) {
      g = &upSymbol;
      *i += 5;
    } else if(strncmp(&s[*i+1], "downkey", 7) == 0) {
      g = &downSymbol;
      *i += 7;
    } else if(strncmp(&s[*i+1], "leftkey", 7) == 0) {
      g = &leftSymbol;
      *i += 7;
    } else if(strncmp(&s[*i+1], "rightkey", 8) == 0) {
      g = &leftSymbol;
      if (hFlip != nullptr) {
        *hFlip = true;
      } else {
        g = &glyphs[0];
      }
      *i += 8;
    } else if(strncmp(&s[*i+1], "leftarrow", 9) == 0) {
      g = &leftArrow;
      *i += 9;
    } else if(strncmp(&s[*i+1], "rightarrow", 10) == 0) {
      g = &rightArrow;
      *i += 10;
    } else if(strncmp(&s[*i+1], "lock", 4) == 0) {
      g = &lockSymbol;
      *i += 4;
    } else if(strncmp(&s[*i+1], "vflip", 5) == 0) {
      *i += 5;
      g = glyphLookup(s, i);
      if (vFlip != nullptr) {
        *vFlip = true;
      }
    } else if (strncmp(&s[*i+1], "hflip", 5) == 0) {
      *i += 5;
      g = glyphLookup(s, i);
      if (vFlip != nullptr) {
        *vFlip = true;
      }
    } else if (strncmp(&s[*i+1], "flat", 4) == 0) {
      g = &flatSymbol;
      *i += 4;
    } else if (strncmp(&s[*i+1], "cent", 4) == 0) {
      g = &centSymbol;
      *i += 4;
    }else {
      g = &glyphs[int(s['\\'])];
    }
  } else {
    g = &glyphs[int(s[*i])];
  }

  return g;
}

/* return number of pixels wide the string is */
int drawChars(ILI9341_t3 *tft, const char *s, uint16_t fgColor, uint16_t midColor, uint16_t bgColor, int justify = justifyLeft) {
  int16_t xOrigin;
  int16_t yOrigin;
  int fullWidth = 0;

  tft->getCursor(&xOrigin, &yOrigin);
  //tft->drawPixel(xOrigin, yOrigin, ILI9341_RED);

  yOrigin -= 5;  /* text is roughly vertically centered on origin */
  
  int i = 0;
  while  (s[i] != '\0') {
    if (s[i] < 0 || s[i] > 127) {
      i++;
      continue;
    }

    CharGlyph* g = glyphLookup(s, &i);
    do {
      fullWidth += g->shift;
      g = g->next;
    } while (g != nullptr);

    fullWidth += charSpacing;
    i++;
  }

  if (fullWidth > 0) {
    fullWidth -= charSpacing;
  }

  if (justify == justifyCenter) {
    xOrigin -= fullWidth/2;
  } else if (justify == justifyRight) {
    xOrigin -= fullWidth;
  }

  i=0;
  while (s[i] != '\0') {
    if (s[i] < 0 || s[i] > 127) {
      continue;
    }

    bool vFlip = false;
    bool hFlip = false;

    CharGlyph *g = glyphLookup(s, &i, &vFlip, &hFlip);
  
    xOrigin += drawGlyph(tft, g, xOrigin, yOrigin, fgColor, midColor, bgColor, vFlip, hFlip);
    xOrigin += charSpacing;

    //Serial.println("draw char " + String(s[i]) + " at " + String(xOrigin) + "," + String(yOrigin));

    i++;
  }

  return fullWidth;
}

#define mkChar(c, width, ...) glyphs[int(c)] = CharGlyph(width, __VA_ARGS__)
#define setShift(c, s) glyphs[int(c)].shift = s

void fontSetup() {
  for (int i=0; i < 128; i++) {
    glyphs[i] = CharGlyph(
      7,
      0b0, 0b1111111,
      0b0, 0b1111111,
      0b0, 0b1111111,
      0b0, 0b1111111,
      0b0, 0b1111111,
      0b0, 0b1111111,
      0b0, 0b1111111,
      0b0, 0b1111111,
      0b0, 0b1111111,
      0b0, 0b1111111);
  }

  mkChar('a', 7,
    0b0000000, 0b0000000,
    0b0000000, 0b0000000,
    0b0000000, 0b0000000,
    0b0000000, 0b0000000,
    0b0011110, 0b0100000,
    0b0100010, 0b1010000,
    0b1000010, 0b0100000, 
    0b1000010, 0b0100000,
    0b0111110, 0b1000001);
  
  setShift('a', 6);
  
  mkChar('b', 6,
    0b100000, 0b000000,
    0b100000, 0b000000,
    0b100000, 0b000000,
    0b100000, 0b000000,
    0b111110, 0b000001,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b111110, 0b000001);

  mkChar('c', 6,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b011111, 0b100000,
    0b100000, 0b000000,
    0b100000, 0b000000,
    0b100000, 0b000000,
    0b011111, 0b100000);

  mkChar('d', 6,
    0b000001, 0b000000,
    0b000001, 0b000000,
    0b000001, 0b000000,
    0b000001, 0b000000,
    0b011111, 0b100000,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b011111, 0b100000);
    
  mkChar('e', 6,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b011110, 0b100001,
    0b100001, 0b000000,
    0b111111, 0b000000,
    0b100000, 0b000000,
    0b011111, 0b100001);

  mkChar('f', 4,
    0b00111, 0b0100,
    0b01000, 0b0000,
    0b01000, 0b0000,
    0b01000, 0b0000,
    0b11110, 0b0000,
    0b01000, 0b0000,
    0b01000, 0b0000,
    0b01000, 0b0000,
    0b01000, 0b0000);

  mkChar('g', 6,
    0b000000, 0b000000,
    0b000000, 0b000000,
    0b000000, 0b000000,
    0b000000, 0b000000,
    0b011111, 0b100000,
    0b100001, 0b010000,
    0b100001, 0b000000,
    0b100001, 0b010000,
    0b011111, 0b100000,
    0b000001, 0b000000,
    0b000001, 0b000000,
    0b111110, 0b000001);

 mkChar('h', 6,
    0b100000, 0b000000,
    0b100000, 0b000000,
    0b100000, 0b000000,
    0b100000, 0b000000,
    0b111110, 0b000001,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b100001, 0b000000);

  mkChar('i', 1,
    0b0, 0b0,
    0b0, 0b0,
    0b1, 0b0,
    0b0, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0);

  mkChar('j', 3,
    0b000, 0b0,
    0b000, 0b0,
    0b001, 0b0,
    0b000, 0b0,
    0b001, 0b0,
    0b001, 0b0,
    0b001, 0b0,
    0b001, 0b0,
    0b001, 0b0,
    0b001, 0b0,
    0b001, 0b0,
    0b110, 0b001);

  mkChar('k', 5,
    0b10000, 0b0,
    0b10000, 0b0,
    0b10000, 0b0,
    0b10000, 0b0,
    0b10011, 0b00100,
    0b11100, 0b0,
    0b10100, 0b00010,
    0b10010, 0b00001,
    0b10001, 0b0);

  mkChar('l', 2,
    0b10, 0b0,
    0b10, 0b0,
    0b10, 0b0,
    0b10, 0b0,
    0b10, 0b0,
    0b10, 0b0,
    0b10, 0b0,
    0b10, 0b0,
    0b10, 0b01);

  setShift('l', 1);

  mkChar('m', 7,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b1110110, 0b0001001,
    0b1001001, 0b0000000,
    0b1001001, 0b0000000,
    0b1001001, 0b0000000,
    0b1001001, 0b0000000);

  mkChar('n', 6,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b111110, 0b000001,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b100001, 0b000000);

  mkChar('o', 6,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b011110, 0b100001,
    0b100001, 0b010010,
    0b100001, 0b000000,
    0b100001, 0b010010,
    0b011110, 0b100001);

  mkChar('p', 6,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b111110, 0b000001,
    0b100001, 0b000010,
    0b100001, 0b000000,
    0b100001, 0b000010,
    0b111110, 0b000001,
    0b100000, 0b0,
    0b100000, 0b0,
    0b100000, 0b0);

  mkChar('q', 8,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b01111100, 0b10000000,
    0b10000100, 0b01000000,
    0b10000100, 0b00000000,
    0b10000100, 0b01000000,
    0b01111100, 0b10000000,
    0b00000100, 0b0,
    0b00000100, 0b00000010,
    0b00000110, 0b00000001);

  setShift('q', 6);

  mkChar('r', 4,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b1011, 0b0100,
    0b1100, 0b0001,
    0b1000, 0b0,
    0b1000, 0b0,
    0b1000, 0b0);

  mkChar('s', 6,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b011111, 0b100000,
    0b100000, 0b000000,
    0b011110, 0b100001,
    0b000001, 0b000000,
    0b111110, 0b100001);

  mkChar('t', 5,
    0b00100, 0b0,
    0b00100, 0b0,
    0b00100, 0b0,
    0b00100, 0b0,
    0b11111, 0b0,
    0b00100, 0b0,
    0b00100, 0b0,
    0b00100, 0b0,
    0b00100, 0b00010);

  mkChar('u', 6,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b100001, 0b0,
    0b100001, 0b0,
    0b100001, 0b0,
    0b100001, 0b010010,
    0b011110, 0b100001);

  mkChar('v', 6,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b100001, 0b0,
    0b100001, 0b0,
    0b100001, 0b010010,
    0b010010, 0b101101,
    0b001100, 0b010010);

  mkChar('w', 7,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b1001001, 0b0,
    0b1001001, 0b0,
    0b1001001, 0b0,
    0b1001001, 0b0110110,
    0b0110110, 0b1001001);

  mkChar('x', 5,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b10001, 0b01010,
    0b01010, 0b10101,
    0b00100, 0b01010,
    0b01010, 0b10101,
    0b10001, 0b01010);

  mkChar('y', 6,
    0b000000, 0b000000,
    0b000000, 0b000000,
    0b000000, 0b000000,
    0b000000, 0b000000,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b100001, 0b010010,
    0b011111, 0b100000,
    0b000001, 0b000000,
    0b000001, 0b000000,
    0b111110, 0b000001);

  mkChar('z', 5,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b11111, 0b0,
    0b00010, 0b00101,
    0b00100, 0b01010,
    0b01000, 0b10100,
    0b11111, 0b0);

mkChar('A', 7,
    0b0111110, 0b1000001,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1111111, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000);

  mkChar('B', 7,
    0b1111100, 0b1000010,
    0b1000010, 0b0000101,
    0b1000001, 0b0000010,
    0b1000001, 0b0000010,
    0b1000010, 0b0000101,
    0b1111110, 0b0000000,
    0b1000001, 0b0000010,
    0b1000001, 0b0000010,
    0b1111110, 0b0000001);

  mkChar('C', 7,
    0b0111110, 0b1000001,
    0b1000001, 0b0100010,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000,
    0b1000001, 0b0100010,
    0b0111110, 0b1000001);

  mkChar('D', 7,
    0b1111100, 0b0000010,
    0b1000010, 0b0000101,
    0b1000001, 0b0000010,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000010,
    0b1000010, 0b0000101,
    0b1111100, 0b0000010);

  mkChar('E', 7,
    0b1111111, 0b1000000,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000,
    0b1111111, 0b0000000,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000,
    0b1111111, 0b0000000);
  
  mkChar('F', 7,
    0b1111111, 0b1000000,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000,
    0b1111111, 0b0000000,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000);

  mkChar('G', 7,
    0b0111110, 0b1000001,
    0b1000001, 0b0100010,
    0b1000000, 0b0000001,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000,
    0b1000111, 0b0001000,
    0b1000001, 0b0000000,
    0b1000001, 0b0100010,
    0b0111110, 0b1000001);

  mkChar('H', 7,
    0b1000001, 0b1000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1111111, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000);

  mkChar('I', 1,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0);

  mkChar('J', 6,
    0b111111, 0b0,
    0b000100, 0b0,
    0b000100, 0b0,
    0b000100, 0b0,
    0b000100, 0b0,
    0b000100, 0b0,
    0b000100, 0b0,
    0b000100, 0b101000,
    0b111000, 0b000100);

  mkChar('K', 7,
    0b1000001, 0b0000010,
    0b1000010, 0b0000101,
    0b1000100, 0b0001010,
    0b1001000, 0b0010100,
    0b1010000, 0b0101000,
    0b1101000, 0b0010100,
    0b1000100, 0b0101010,
    0b1000010, 0b0000101,
    0b1000001, 0b0000010);

  mkChar ('L', 6,
    0b100000, 0b0,
    0b100000, 0b0,
    0b100000, 0b0,
    0b100000, 0b0,
    0b100000, 0b0,
    0b100000, 0b0,
    0b100000, 0b0,
    0b100000, 0b000001,
    0b111111, 0b0);

  mkChar('M', 7,
    0b1000001, 0b0100010,
    0b1100011, 0b00,
    0b1100011, 0b0010100,
    0b1010101, 0b0100010,
    0b1010101, 0b0001000,
    0b1001001, 0b0,
    0b1000001, 0b0001000,
    0b1000001, 0b0,
    0b1000001, 0b0);

  mkChar('N', 6,
    0b100001, 0b010000,
    0b110001, 0b001001,
    0b101001, 0b010100,
    0b100101, 0b001010,
    0b100011, 0b000100,
    0b100001, 0b000010,
    0b100001, 0b0,
    0b100001, 0b0,
    0b100001, 0b0);

  mkChar('O', 7,
    0b0111110, 0b1000001,
    0b1000001, 0b0100010,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0100010,
    0b0111110, 0b1000001);

  mkChar('P', 7,
    0b1111110, 0b1000001,
    0b1000001, 0b0000010,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000010,
    0b1111110, 0b0000001,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000,
    0b1000000, 0b0000000);

  mkChar('Q', 8,
    0b01111100, 0b10000010,
    0b10000010, 0b01000100,
    0b10000010, 0b00000000,
    0b10000010, 0b00000000,
    0b10000010, 0b00000000,
    0b10000010, 0b00000000,
    0b10011010, 0b00000100,
    0b10000110, 0b01000000,
    0b01111101, 0b10000010);

  mkChar('R', 7,
    0b1111110, 0b1000001,
    0b1000001, 0b0000010,
    0b1000001, 0b0000000,
    0b1000001, 0b0000000,
    0b1000001, 0b0000010,
    0b1111110, 0b0000001,
    0b1000010, 0b0000100,
    0b1000001, 0b0000010,
    0b1000001, 0b0000000);

  mkChar('S', 7,
    0b0111111, 0b1000000,
    0b1000000, 0b0100000,
    0b1000000, 0b0000000,
    0b1000000, 0b0100000,
    0b0111110, 0b1000001,
    0b0000001, 0b0000010,
    0b0000001, 0b0000000,
    0b0000001, 0b0000010,
    0b1111110, 0b0000001);

  mkChar('T', 7,
    0b1111111, 0b0,
    0b0001000, 0b0,
    0b0001000, 0b0,
    0b0001000, 0b0,
    0b0001000, 0b0,
    0b0001000, 0b0,
    0b0001000, 0b0,
    0b0001000, 0b0,
    0b0001000, 0b00101000);

  mkChar('U', 6,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b100001, 0b000000,
    0b100001, 0b010010,
    0b011110, 0b100001);

  mkChar('V', 6,
    0b100001, 0b0000000,
    0b100001, 0b0000000,
    0b100001, 0b0000000,
    0b100001, 0b0000000,
    0b100001, 0b0000000,
    0b100001, 0b0100010,
    0b010010, 0b1000001,
    0b010010, 0b0011000,
    0b001100, 0b0000000);

  /* w is 9 pixels wide, so we have to represent it as a double glyph */
  mkChar('W', 8,
    0b10000000, 0b0,
    0b10000000, 0b0,
    0b10000000, 0b0,
    0b10000000, 0b0,
    0b01000001, 0b10010000,
    0b01001001, 0b10101000,
    0b01010101, 0b00010000,
    0b00110110, 0b01000001,
    0b00100010, 0b00000000);

  wExtra = CharGlyph(1,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b0, 0b1,
    0b0, 0b1,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0);

  glyphs[int('W')].next = &wExtra;

  mkChar('X', 7,
    0b1000001, 0b0000000,
    0b1000001, 0b0100010,
    0b0100010, 0b1010101,
    0b0010100, 0b0101010,
    0b0001000, 0b0010100,
    0b0010100, 0b0101010,
    0b0100010, 0b1010101,
    0b1000001, 0b0100010,
    0b1000001, 0b0000000);

  mkChar('Y', 7,
    0b1000001, 0b0,
    0b1000001, 0b0100010,
    0b0100010, 0b0010100,
    0b0010100, 0b0101010,
    0b0001000, 0b0010100,
    0b0001000, 0b0,
    0b0001000, 0b0,
    0b0001000, 0b0,
    0b0001000, 0b0);

 mkChar('Z', 7,
    0b1111111, 0b0000000,
    0b0000001, 0b0000010,
    0b0000010, 0b0000101,
    0b0000100, 0b0001010,
    0b0001000, 0b0010100,
    0b0111110, 0b1000001,
    0b0100000, 0b1010000,
    0b1000000, 0b0100000,
    0b1111111, 0b0000000);

  mkChar(' ', 5,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0);

  mkChar('.', 1,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b1, 0b0);

  mkChar(':', 1,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b1, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b1, 0b0,
    0b0, 0b0,
    0b0, 0b0);

  mkChar('-', 4,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0110, 0b1001,
    0b0, 0b0,
    0b0, 0b0);

  mkChar('@', 8,
    0b00000000, 0b00000000,
    0b00000000, 0b00000000,
    0b01111111, 0b10000000,
    0b10000000, 0b01000001,
    0b10001110, 0b00010000,
    0b10010010, 0b00101000,
    0b10100010, 0b00010000, 
    0b10100010, 0b00000000,
    0b10011111, 0b00100000);
  
  atExtra = CharGlyph(1,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b1,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0);

  glyphs[int('@')].next = &atExtra;

  mkChar('>', 4,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b10000, 0b01000,
    0b01000, 0b10100,
    0b00100, 0b01010,
    0b00010, 0b00101,
    0b00100, 0b01010,
    0b01000, 0b10100,
    0b10000, 0b01000);

  mkChar('<', 5,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b00001, 0b00010,
    0b00010, 0b00101,
    0b00100, 0b01010,
    0b01000, 0b10100,
    0b00100, 0b01010,
    0b00010, 0b00101,
    0b00001, 0b00010);

  mkChar('(', 3,

    0b001, 0b010,
    0b010, 0b001,
    0b010, 0b100,
    0b100, 0b010,
    0b100, 0b000,
    0b100, 0b000,
    0b100, 0b010,
    0b010, 0b100,
    0b010, 0b001,
    0b001, 0b010);

  mkChar(')', 3,
   
    0b100, 0b010,
    0b010, 0b100,
    0b010, 0b001,
    0b001, 0b010,
    0b001, 0b000,
    0b001, 0b000,
    0b001, 0b010,
    0b010, 0b001,
    0b010, 0b100,
    0b100, 0b010);
  
  mkChar('+', 5,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b00000, 0b00100,
    0b00100, 0b00000,
    0b01110, 0b10001,
    0b00100, 0b00000,
    0b00000, 0b00100);

  mkChar('=', 5,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b00000, 0b00000,
    0b01110, 0b10001,
    0b00000, 0b00000,
    0b01110, 0b10001,
    0b00000, 0b00000);

  mkChar('*', 5,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b00100, 0b00000,
    0b10101, 0b01010,
    0b01110, 0b00000,
    0b10101, 0b01010,
    0b00100, 0b00000,
    0b0, 0b0);

  mkChar('_', 6,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b011110, 0b100001);

  mkChar('\'', 1,
    0b0, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b0, 0b1,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0);

  mkChar('\"', 3,
    0b0, 0b0,
    0b101, 0b000,
    0b101, 0b000,
    0b000, 0b101,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0);

  mkChar('~',7,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0,
    0b0000000, 0b0110000,
    0b0111110, 0b1000001,
    0b0000000, 0b0000110,
    0b0, 0b0,
    0b0, 0b0,
    0b0, 0b0);

  mkChar('0', 6,
    0b011110, 0b100001,
    0b100001, 0b010010,
    0b100011, 0b000100,
    0b100101, 0b001010,
    0b101001, 0b010100,
    0b110001, 0b001000,
    0b100001, 0b010000,
    0b100001, 0b010010,
    0b011110, 0b100001);

  mkChar('1', 1,
    0b01, 0b0,
    0b01, 0b10,
    0b01, 0b0,
    0b01, 0b0,
    0b01, 0b0,
    0b01, 0b0,
    0b01, 0b0,
    0b01, 0b0,
    0b01, 0b0);

  mkChar('2', 6,
    0b011110, 0b100001,
    0b100001, 0b010010,
    0b000001, 0b100010,
    0b000010, 0b000101,
    0b000100, 0b001010,
    0b011000, 0b100100,
    0b100000, 0b010000,
    0b100000, 0b0,
    0b111111, 0b0);

  mkChar('3', 5,
    0b11110, 0b00001,
    0b00001, 0b00010,
    0b00001, 0b00000,
    0b00001, 0b00010,
    0b01110, 0b00001,
    0b00001, 0b00010,
    0b00001, 0b00000,
    0b00001, 0b00000,
    0b11110, 0b00001);

  mkChar('4', 6,
    0b000010, 0b100000,
    0b100010, 0b0,
    0b100010, 0b0,
    0b100010, 0b0,
    0b100010, 0b000001,
    0b111110, 0b000001,
    0b000010, 0b0,
    0b000010, 0b0,
    0b000010, 0b0);

  setShift('4', 5);

  mkChar('5', 6,
    0b111111, 0b000000,
    0b100000, 0b000000,
    0b100000, 0b000000,
    0b111110, 0b000001,
    0b000001, 0b000010,
    0b000001, 0b000000,
    0b000001, 0b000000,
    0b000001, 0b000010,
    0b011110, 0b100001);

  mkChar('6', 6,
    0b011111, 0b100000,
    0b100000, 0b010001,
    0b100000, 0b000000,
    0b100000, 0b010000,
    0b011110, 0b100001,
    0b100001, 0b010010,
    0b100001, 0b000000,
    0b100001, 0b010010,
    0b011110, 0b100001);

  mkChar('7', 5,
    0b11111, 0b0,
    0b00001, 0b10010,
    0b00010, 0b00001,
    0b00010, 0b00100,
    0b11111, 0b00000,
    0b00100, 0b01000,
    0b01000, 0b00100,
    0b01000, 0b00000,
    0b01000, 0b00000);

  mkChar('8', 6,
    0b011110, 0b100001,
    0b100001, 0b010010,
    0b100001, 0b000000,
    0b100001, 0b010010,
    0b011110, 0b100001,
    0b100001, 0b010010,
    0b100001, 0b000000,
    0b100001, 0b010010,
    0b011110, 0b100001);

  mkChar('9', 6,
    0b011110, 0b100001,
    0b100001, 0b010010,
    0b100001, 0b000000,
    0b100001, 0b010010,
    0b011110, 0b100001,
    0b000001, 0b000010,
    0b000001, 0b000000,
    0b000001, 0b100010,
    0b111110, 0b000001);

  mkChar('#', 8,
    0b00000100, 0b00100000,
    0b00100111, 0b00001000,
    0b00111100, 0b01000010,
    0b11100100, 0b00010000,
    0b00100111, 0b00001000,
    0b00111100, 0b01000010,
    0b11100100, 0b00010000,
    0b00100100, 0b00000000,
    0b00100000, 0b00000100);

  mkChar('!', 1,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b1, 0b0,
    0b0, 0b0,
    0b1, 0b0);

  mkChar('/', 4,
    0b0000, 0b0000,
    0b0001, 0b0000,
    0b0001, 0b0010,
    0b0010, 0b0001,
    0b0010, 0b0100,
    0b0100, 0b0010,
    0b0100, 0b1000,
    0b1000, 0b0100,
    0b1000, 0b0000);

  mkChar('%', 6,
    0b000000, 0b000000,
    0b000010, 0b000000,
    0b100010, 0b000100,
    0b000100, 0b000010,
    0b000100, 0b001000,
    0b001000, 0b000100,
    0b001000, 0b010000,
    0b010001, 0b001000,
    0b010000, 0b000000);

  din5Plug.next = &din5Plug2;
  din5PlugSmall.next = &din5PlugSmall2;

  usbSymbol.next = &usbSymbol2;
  usbSymbol2.next = &usbSymbol3;
  usbSymbol3.next = &usbSymbol4;

  upSymbol.next = &upSymbol2;
  upSymbol.flipUpDown();
  upSymbol2.flipUpDown();

  downSymbol.next = &downSymbol2;
  leftSymbol.next = &leftSymbol2;
}
