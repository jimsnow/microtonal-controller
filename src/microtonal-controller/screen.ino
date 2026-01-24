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

/* Screen */

#define TFT_DC      screenDCPin
#define TFT_CS      screenCSPin
#define TFT_RST     255  // 255 = unused, connect to 3.3V
#define TFT_MOSI    sdiPin
#define TFT_SCLK    sckPin
#define TFT_MISO    sdoPin
ILI9341_t3 tft = ILI9341_t3(TFT_CS, TFT_DC, TFT_RST, TFT_MOSI, TFT_SCLK, TFT_MISO);

uint16_t darkGrey = tft.color565(48, 48, 48);

//#define width 320
const int width = 320;
#define menuSpacing 2
#define height 240
#define menuWidth 130
#define menuItemHeight (48-(4*menuSpacing))
#define menuItemSpacing (5*menuSpacing)
#define statusWidth (width - (menuWidth + 16))
#define statusHeight 20
#define navButtonWidth ((statusWidth / 2) - 6)
#define visualizerWidth (statusWidth)
#define visualizerHeight (height - (menuItemHeight * 2 + statusHeight * 2) + menuItemSpacing)

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
  String textL;
  String textR;
  String textC;
  uint16_t bgcolor;
  uint16_t fgcolor;
  bool enabled;
  bool redraw;
  bool highlight;
  uint32_t textUsecs;
};

struct Point cursorOld(const struct Window &window) {
  uint16_t left = window.extent.p1.x;
  //uint16_t right = window.extent.p2.x;
  uint16_t top = window.extent.p1.y;
  uint16_t bottom = window.extent.p2.y;
  uint16_t bheight = bottom-top;

  return Point(left+5, top+(bheight/2) - 9);
}

struct Point cursor(const struct Window &window, int justify) {
  uint16_t left = window.extent.p1.x;
  uint16_t right = window.extent.p2.x;
  uint16_t top = window.extent.p1.y;
  uint16_t bottom = window.extent.p2.y;
  uint16_t y = (top + bottom) / 2;

  switch (justify) {
    case(justifyLeft):   return Point(left+10, y);
    case(justifyCenter): return Point((left+right) / 2, y);
    case(justifyRight):  return Point(right-10, y);
  }

  return Point(left+5, (top + bottom) / 2);
}

void setWindowCursor(const struct Window &window, int justify) {
  Point c = cursor(window, justify);
  tft.setCursor(c.x, c.y);
}

bool moreUp = false;
bool moreDown = false;

#define numWindows 12
struct Window windows[numWindows];

void renderScreen(uint32_t deltaUsecs) {
  for (int i=0; i<numWindows; i++) {
    if (!windows[i].redraw) {
      continue;
    }

    Rectangle r = windows[i].extent;
    if (!windows[i].enabled) {
      tft.fillRect(r.p1.x, r.p1.y, (r.p2.x-r.p1.x), (r.p2.y-r.p1.y), ILI9341_BLACK);
      tft.fillRect(r.p1.x, r.p1.y+2, (r.p2.x-r.p1.x), (r.p2.y-r.p1.y)-4, windows[i].bgcolor);
      tft.fillRect(r.p1.x+1, r.p1.y+3, (r.p2.x-r.p1.x)-2, (r.p2.y-r.p1.y)-6, ILI9341_BLACK);
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

    tft.fillRect(r.p1.x, r.p1.y+2, (r.p2.x-r.p1.x), (r.p2.y-r.p1.y)-2, bgcolor);
    
    /*
    tft.setTextColor(fgcolor);
    tft.setTextSize(2);
    setWindowCursor(windows[i]);
    tft.println(windows[i].text);
    */

    setWindowCursor(windows[i], justifyLeft);
    drawChars(&tft, windows[i].textL.c_str(), fgcolor, ILI9341_LIGHTGREY, bgcolor, justifyLeft);

    setWindowCursor(windows[i], justifyCenter);
    drawChars(&tft, windows[i].textC.c_str(), fgcolor, ILI9341_LIGHTGREY, bgcolor, justifyCenter);

    setWindowCursor(windows[i], justifyRight);
    drawChars(&tft, windows[i].textR.c_str(), fgcolor, ILI9341_LIGHTGREY, bgcolor, justifyRight);

    windows[i].redraw = false;
  }

  if (moreUp) {
    //tft.setTextColor(windows[menuText1].highlight ? windows[menuText1].bgcolor : windows[menuText1].fgcolor);
    tft.setCursor(menuWidth/2, 10);
    //tft.println("...");

    if (windows[menuText1].highlight) {
      drawChars(&tft, "\\upkey", windows[menuText1].bgcolor, ILI9341_LIGHTGREY, windows[menuText1].fgcolor, justifyCenter);
    } else {
      drawChars(&tft, "\\upkey", windows[menuText1].fgcolor, ILI9341_LIGHTGREY, windows[menuText1].bgcolor, justifyCenter);
    }
  }

  if (moreDown) {
    //tft.setTextColor(windows[menuText5].highlight ? windows[menuText5].bgcolor : windows[menuText5].fgcolor);
    tft.setCursor(menuWidth/2, height-8);
    //tft.println("...");

    if (windows[menuText5].highlight) {
      drawChars(&tft, "\\downkey", windows[menuText1].bgcolor, ILI9341_LIGHTGREY, windows[menuText1].fgcolor, justifyCenter);
    } else {
      drawChars(&tft, "\\downkey", windows[menuText1].fgcolor, ILI9341_LIGHTGREY, windows[menuText1].bgcolor, justifyCenter);
    }

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
    windows[i].textL = "";
    windows[i].textC = "";
    windows[i].textR = "";
    windows[i].bgcolor = darkGrey;
    windows[i].fgcolor = ILI9341_WHITE;
    windows[i].enabled = true;
    windows[i].redraw = true;
    windows[i].highlight = false;
    windows[i].textUsecs = 0;
  }


  for (int i=0; i<5; i++) {
    windows[i].extent = Rectangle (0, i*menuItemHeight + i*menuItemSpacing, menuWidth, (i+1)*menuItemHeight + i*menuItemSpacing);
  }

  windows[backText].extent  = Rectangle ((width - statusWidth), 0, (width - statusWidth) + navButtonWidth, menuItemHeight);
  windows[fwdText].extent  = Rectangle (width - navButtonWidth, 0, width, menuItemHeight);
  windows[cancelText].extent  = Rectangle ((width - statusWidth), menuItemHeight + menuItemSpacing, (width - statusWidth) + navButtonWidth, menuItemHeight * 2 + menuItemSpacing);
  windows[okText].extent  = Rectangle (width - navButtonWidth, menuItemHeight + menuItemSpacing, width, menuItemHeight * 2 + menuItemSpacing);

  windows[statusBar1].extent = Rectangle (width-statusWidth, menuItemHeight * 2 + menuItemSpacing * 2, width, menuItemHeight * 2 + menuItemSpacing * 2 + statusHeight);
  windows[statusBar1].bgcolor = ILI9341_BLACK;
  windows[statusBar2].extent  = Rectangle (width-statusWidth, height-statusHeight, width, height);
  windows[statusBar2].bgcolor = ILI9341_BLACK;

  windows[visualizerWindow].extent = Rectangle (width-statusWidth, menuItemHeight * 2 + statusHeight, width, height - statusHeight);
  windows[visualizerWindow].bgcolor = ILI9341_BLACK;

  windows[backText].textC = "back";
  windows[backText].enabled = false;
  windows[fwdText].enabled = false;
  windows[cancelText].enabled = false;
  windows[okText].enabled = false;

  pinMode(backlightPin, OUTPUT);
  analogWriteFrequency(backlightPin, 3611*2); /* default is 3.611 kHz */
  analogWrite(backlightPin, brightness);
  brightnessSet = brightness;

  delayMicroseconds(100000);

  tft.begin();
  delayMicroseconds(100000);
  tft.setRotation(1);
  tft.setClock(100000000);
  tft.fillScreen(ILI9341_BLACK);

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
  /* if something else is using the visualizerWindow, don't update */
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

/* Menu */

/* menu item that's selected for editing */
struct MenuItem* editItem = nullptr;


#define menuStackSize 10
struct MenuItem* menuStack[menuStackSize];
uint16_t menuStackPos = 0; /* points to first unoccupied slot */


struct MenuItem* menu[9] = {&emptyMenuItem, &emptyMenuItem, &emptyMenuItem, &emptyMenuItem, &emptyMenuItem, &emptyMenuItem, &emptyMenuItem,  &emptyMenuItem, &emptyMenuItem};


void menuSelect(struct MenuItem *item, uint16_t button) {
  if (item == nullptr) {
    return;
  }

  if (item->type == submenu) {
    for (int i=0; i < 5; i++) {
      if (i + item->scrollOffset < item->numChildren) {
        //windows[i].bgcolor = ILI9341_DARKGREY;
        windows[i].textL = item->children[i]->text;
        windows[i].enabled = item->children[i]->type != empty;
        menu[i] = item->children[i];
        Serial.println("displaying menu item " + String(i) + " " + windows[i].textL);
      } else {
        windows[i].enabled = false;
        menu[i] = &emptyMenuItem;
      }
      windows[i].redraw = true;
    }

    for (int i=0; i < 3; i++) {
      if (item->options[i] != nullptr) {
        //windows[6+i].bgcolor = ILI9341_DARKGREY;
        windows[6+i].textC = item->options[i]->text;
        windows[6+i].enabled = item->options[i]->type != empty;
        menu[6+i] = item->options[i];
      } else {
        windows[6+i].enabled = false;
        menu[6+i] = &emptyMenuItem;
      }
      windows[6+i].redraw = true;
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

  if (item->type == toggle) {
    if (item->data != nullptr) {
      *((bool*)(item->data)) = !*(bool*)(item->data);
      Serial.println("toggled " + String(*((bool*)(item->data))));
    }
  }

  if (item->type == selection) {
    if (item->data != nullptr) {
      *((uint32_t*)(item->data)) = item->defaultValue;
    }
  }

  if (item->type == value) {
    windows[visualizerWindow].textC = "\\leftkey " + String(*(uint32_t *)(item->data)) + " \\rightkey";
    windows[visualizerWindow].enabled = true;
    editItem = item;
    windows[visualizerWindow].redraw = true;
  } else if (item->type == floatValue) {
    windows[visualizerWindow].textC = "\\leftkey " + String(*(float*)(item->data), 3) + " \\rightkey";
    windows[visualizerWindow].enabled = true;
    editItem = item;
    windows[visualizerWindow].redraw = true;
  } else if (windows[visualizerWindow].enabled == true) {
    editItem = nullptr;
    windows[visualizerWindow].textL = "";
    windows[visualizerWindow].textC = "";
    windows[visualizerWindow].textR = "";
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

      windows[i].highlight = child->checkHighlight(editItem);
      if (windows[i].highlight != prevHighlight) {
        windows[i].redraw = true;
      }
    }

    for (int i = 0; i < 3; i++) {
      auto option = parent->options[i];
      if (option == nullptr) {
        continue;
      }
      bool prevHighlight = option->highlight;
      windows[6+i].highlight = option->checkHighlight(editItem);
      if (windows[6+i].highlight != prevHighlight) {
        windows[6+i].redraw = true;
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
  if (menuStackPos == 0) {
    return;
  }

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
      windows[i].textL = menu[i]->text;
      windows[i].highlight = item->children[i]->checkHighlight(editItem);
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
  if (button == backText) {
    Serial.println("back button released");
    menuBack();
  } else {
    Serial.println("released button " + String(button));
    menuSelect(menu[button], button);
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
      status1->textL = "";
      status1->textC = "";
      status1->textR = "";
      status1->redraw = true;
    } else {
      status1->textUsecs -= deltaUsecs;
    }
  }
}

void status1TextUpdate(String textC, uint32_t usecs) {
  windows[statusBar1].textC = textC;
  windows[statusBar1].textUsecs = usecs;
  windows[statusBar1].redraw = true;
}

void setStatusText(String textL, String textC, String textR) {
  windows[statusBar2].textL = textL;
  windows[statusBar2].textC = textC;
  windows[statusBar2].textR = textR;
  windows[statusBar2].redraw = true;
}

void incrementMenuValue() {
  if (editItem == nullptr) {
    return;
  }

  if (editItem->type == value) {
    incrementIntValue();
  } else if (editItem->type == floatValue) {
    incrementFloatValue();
  }
}

void decrementMenuValue() {
  if (editItem == nullptr) {
    return;
  }

  if (editItem->type == value) {
    decrementIntValue();
  } else if (editItem->type == floatValue) {
    decrementFloatValue();
  }
}

void incrementIntValue() {
  auto val = *(uint32_t*)(editItem->data);

  if (editItem->maxValue == nullptr || val < *editItem->maxValue) {
    val++;
    windows[visualizerWindow].textC = "\\leftkey " + String(val) + " \\rightkey";
    *(uint32_t*)(editItem->data) = val;
    windows[visualizerWindow].redraw = true;
  }
}

void decrementIntValue() {
  auto val = *(uint32_t*)(editItem->data);

  if (editItem->minValue == nullptr || val > *editItem->minValue) {
    val--;
    windows[visualizerWindow].textC = "\\leftkey " + String(val) + " \\rightkey";
    *(uint32_t*)(editItem->data) = val;
    windows[visualizerWindow].redraw = true;
  }
}

void incrementFloatValue() {
  auto val = *(float*)(editItem->data);

  if (val < 1.0f) {
    val += 0.005;
    val = clamp(val);
    windows[visualizerWindow].textC = "\\leftkey " + String(val, 3) + " \\rightkey";
    *(float*)(editItem->data) = val;
    windows[visualizerWindow].redraw = true;
  }
}

void decrementFloatValue() {
  auto val = *(float*)(editItem->data);

  if (val > 0.0f) {
    val -= 0.005;
    val = clamp(val);
    windows[visualizerWindow].textC = "\\leftkey " + String(val, 3) + " \\rightkey";
    *(float*)(editItem->data) = val;
    windows[visualizerWindow].redraw = true;
  }
}

void redrawVisualizer() {
  windows[visualizerWindow].redraw = true;
}

void fontTest() {
  windows[visualizerWindow].enabled = true;

  setWindowCursor(windows[visualizerWindow], justifyCenter);

  int16_t xOrigin, yOrigin;
  tft.getCursor(&xOrigin, &yOrigin);

  for (int i = 2; i < 8; i++) {
    for (int j = 0; j < 16; j++) {
      int c = i * 16 + j;
      int x = xOrigin + (j - 8) * 10 + 5;
      int y = yOrigin + (i - 5) * 14;

      drawGlyph(&tft, &glyphs[c], x, y, windows[visualizerWindow].fgcolor, ILI9341_LIGHTGREY, windows[visualizerWindow].bgcolor);
    }
  }

  // don't set redraw here...
}
