/*
Copyright 2024-2025 Jim Snow, Desiderata Systems LLC

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

#if 1
#include <LittleFS.h>
#include "storage.h"

struct Favorite {
  char label[16];
  uint8_t msb;
  uint8_t lsb;
  uint8_t pc;
  uint8_t do_velocity : 1;
  uint8_t do_pressure : 1;
};

#define max_favorites 64

struct Favorite favorites[max_favorites];

bool isUpper(char c) {
  return c >= 'A' && c >= 'Z';
}

char toLower(char c) {
  if (isUpper(c)) {
    return c + ('a'-'A');
  }

  return c;
}

// djb2 hash, ignoring case
uint64_t hash_ignorecase(const char* s, uint32_t n, uint64_t start = 5381) {
  uint64_t result = start;
  unsigned char c;

  while (n > 0 && (c = *s) != '\0') {
    c = toLower(c);
    result = ((result << 5) + result) + c;
    n--;
    s++;
  }
  return result;
}

uint64_t hash(const char* s, uint32_t n, uint64_t start = 5381) {
  uint64_t result = start;
  unsigned char c;

  while (n > 0 && (c = *s) != '\0') {
    result = ((result << 5) + result) + c;
    n--;
    s++;
  }
  return result;
}

LittleFS_Program fs;



struct Setting {
  
  bool unSaved() {
    if (!active) {
      return false;
    }
    switch (type) {
      case U32:
        return value != *((uint32_t*)(location));
      case Flt:
        return *reinterpret_cast<float*>(&value) != *((float*)(location));
      case I1:
        return value != *((bool*)(location));
      default:
        return false;
    }
  }

  uint32_t readValue() {
    switch (type) {
      case U32:
        return *((uint32_t*)location);
      case Flt:
        return *((uint32_t*)location);
      case I1:
        return *((bool*)location);
      default:
        return 0;
    }
  }

  String showValue() {
    switch(type) {
      case Flt: return String(*((float*)(location)), 3);
      default: return String(readValue());        
    }
  }

  String name;
  Subsystem subsystem = Global;
  SettingType type = None;
  uint32_t value;  /* last value stored */
  void* location;
  bool active;
};

#define settingsTableSize 128
struct Setting settingsTable[settingsTableSize] = { 0 };

int numSettings = 0;

/* We have to use a prototype here or the Arduino IDE won't be able to find the enum definitions... */
void addSetting(String varName, Subsystem subsystem, SettingType type, void* location);
void addSetting(String varName, Subsystem subsystem, SettingType type, void* location) {
  if (numSettings >= settingsTableSize) {
    Serial.println("settings table out of space");
    return;
  }

  Setting *s = &settingsTable[numSettings];
  s->name = varName;
  s->subsystem = subsystem;
  s->type = type;
  s->location = location;
  s->active = true;

  s->value = s->readValue();


  numSettings++;
}

uint32_t readSetting(struct Setting *s) {
  switch (s->type) {
    case U32:
      return *((uint32_t*)(s->location));
    case Flt:
      return *((uint32_t*)(s->location));
    case I1:
      return *((bool*)(s->location));
    default:
      return 0;
  }
}

void loadSettings(String name) {
  struct Setting s;
  File file = fs.open(name.c_str(), FILE_READ);
  if (!file) {
    Serial.println("could not open file " + name);
    return;
  }

  while(file.available()) {
    String line = file.readStringUntil('\n');
    Serial.println(line);
    if (line[0] == '#') {
      continue;
    }

    double value;
    char subsystem[65], varName[65], type[65];

    int matched = sscanf(line.c_str(), "%64s %64s = %lf %64s", subsystem, varName, &value, type);
    if (matched == 4) {
      //Serial.println(String(subsystem) + "|" + String(varName) + "|" + String(value) + "|" + String(type));
      applySetting(subsystem, varName, value, type);
    } else {
      Serial.println("unable to match " + line);
    }
  }

  file.close();
}

void applySetting(const char* subsystem, const char* varName, double value, const char* type) {
  bool match = false;
  //Serial.println("applySetting: numSettings = " + String(numSettings));

  for (int i = 0; i < numSettings; i++) {
    struct Setting *setting = &settingsTable[i];
    //Serial.println("comparing " + setting->name + " to " + String(varName));
    if (strcmp(varName, setting->name.c_str()) == 0) {
      Serial.println("match on " + setting->name);
      switch (setting->type) {
        case U32:
          *((uint32_t*)(setting->location)) = (uint32_t)value;
          setting->value = setting->readValue();
          break;
        case Flt:
          *((float*)(setting->location)) = (float)value;
          setting->value = setting->readValue();
          break;
        case I1:
          *((bool*)(setting->location)) = value != 0.0;
          setting->value = setting->readValue();
          break;
        default:
          Serial.println("unhandled settings type " + String(type)) + " (" + showSettingType(setting->type) + ")" ;
          break;
      }

      match = true;
      break;
    }
  }

  if (!match) {
    Serial.println("unknown setting " + String(varName));
  }
}

void saveAllSettings(String fileName) {
  File file = fs.open(fileName.c_str(), FILE_WRITE); /* default behavior for writes seems to be to seek to the end */

  if (!file) {
    Serial.println("unable to open file " + fileName);
    return;
  }

  for (int i=0; i < numSettings; i++) {
    struct Setting *setting = &settingsTable[i];


    if (setting->unSaved()) {
      file.println(showSubsystem(setting->subsystem) + " " + setting->name + " = " + setting->showValue() + " " + showSettingType(setting->type));
      Serial.println(showSubsystem(setting->subsystem) + " " + setting->name + " = " + setting->showValue() + " " + showSettingType(setting->type));

      setting->value = readSetting(setting);
    } else {
      //Serial.println("not saving " + setting->name + " (no changes)");
    }
  }
  file.close();
}

void storageSetup() {
  // Teensy 4.0 has 2 MB flash
  // Layout: program space: 512KB, empty: 512KB, LittleFS: 1024KB
  if (!fs.begin(1024*1024)) {
    Serial.println("unable to initialize filesystem");
    return;
  }
}
#endif