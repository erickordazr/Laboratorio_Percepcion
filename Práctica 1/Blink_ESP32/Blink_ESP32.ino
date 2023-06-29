#define LED_BUILTIN 5

void setup() {
// Inicializar el pin digital LED_BUILTIN como una salida.
pinMode(LED_BUILTIN, OUTPUT);
}

// La función loop se ejecuta una y otra vez para siempre.
void loop() {
digitalWrite(LED_BUILTIN, HIGH); // Encender el LED (HIGH es el nivel de voltaje)
delay(1000); // Esperar un segundo
digitalWrite(LED_BUILTIN, LOW); // Apagar el LED haciendo el voltaje LOW
delay(1000); // Esperar un segundo
}
