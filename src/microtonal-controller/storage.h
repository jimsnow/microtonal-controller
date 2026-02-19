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

#ifndef STORAGE_H
#define STORAGE_H

enum SettingType {
  None = 0,
  I1,
  I8,
  I16,
  I32,
  U8,
  U16,
  U32,
  Flt
};

static inline String showSettingType(SettingType t) {
  switch (t) {
    case None: return "none";
    case I1: return "bool";
    case I32: return "int32_t";
    case U32: return "uint32_t";
    case Flt: return "float";
    default: return "unknown";
  }
}

enum Subsystem {
  Global,
  Controls,
  Midi,
  Synth,
  Tuning,
  Interface
};

static inline String showSubsystem(Subsystem s) {
  switch (s) {
    case Global: return "global";
    case Controls: return "controls";
    case Midi: return "midi";
    case Synth: return "synth";
    case Tuning: return "tuning";
    case Interface: return "interface";
    default: return "invalid";
  }
}

/* We have to use a prototype here or the Arduino IDE won't be able to find the enum definitions... */
void addSetting(String varName, Subsystem subsystem, SettingType type, void* location);
#endif