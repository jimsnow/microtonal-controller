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
  MenuItem(String text, enum menuItemType type, uint32_t *data, uint32_t *minValue, uint32_t *maxValue) : text{text}, type{type}, data{data}, minValue{minValue}, maxValue{maxValue} {}
  MenuItem(String text, float *data) : text{text}, type{floatValue}, data{(void*)data} {}
  MenuItem(String text, bool *data) : text{text}, type{toggle}, data{(void*)data} {
    if (data != nullptr) {
      highlight = *((bool*)data);
    }
  }
  MenuItem(String text, enum menuItemType type, void *data, uint32_t defaultValue = 0, void (*select)(void *data) = nullptr) : text{text}, type{type}, select{select}, data{data}, defaultValue{defaultValue} {
    if (type == selection) {
      highlight = defaultValue == *(uint32_t*)data;
    }
  }


  bool checkHighlight(struct MenuItem *editItem) {
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
      case floatValue:
        highlight = editItem == this;
        break;
      default:
        highlight = false;
        break;
    }

    return highlight;
  }

  void addOption(struct MenuItem *item, int index) {
    if (index < 0 || index > 2) {
      return;
    }

    options[index] = item;
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
  struct MenuItem* options[3] = {};
};
