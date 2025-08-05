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


FlexCAN_T4<CAN1, RX_SIZE_256, TX_SIZE_16> can;

void canSetup() {
  can.begin();
  can.setBaudRate(1000000);
}

void canUpdate() {
  CAN_message_t msg;
  if (can.read(msg)){
    Serial.print("can message received ID 0x");
    Serial.print(msg.id, HEX);
    Serial.print(" data ");
    for (int i = 0; i < 8; i++) {
      Serial.print(msg.buf[i], HEX); Serial.print(" ");
    }
    Serial.println();
  }

  /*
  static uint32_t t_start = millis();
  if (millis() > t_start + 100) {
    Serial.println("sending ping");
    msg.id = random(0x1,0x7FE);
    for ( uint8_t i = 0; i < 8; i++ ) msg.buf[i] = i + 1;
    can.write(msg);
    t_start = millis();
  } */
}
