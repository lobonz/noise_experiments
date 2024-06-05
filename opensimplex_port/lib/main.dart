import 'dart:io';
import 'dart:math';
import 'dart:ui' as ui;
import 'package:flutter/material.dart';
import 'package:opensimplex_port/open_simplex2d.dart';
import 'package:path_provider/path_provider.dart';
import 'package:flutter/services.dart';

class OpenSimplexDemo extends StatelessWidget {
  static const int width = 100;
  static const int height = 100;
  static const double frequency = 1.0 / 24.0;
  static const int seed = 0;

  const OpenSimplexDemo({super.key});

  Future<void> generateNoiseImage() async {
    final ByteData byteData = await rootBundle.load('assets/noise.png');
    final Uint8List uint8List = byteData.buffer.asUint8List();

    ui.Image noiseImage = await decodeImageFromList(uint8List);
    ByteData? pngBytes =
        await noiseImage.toByteData(format: ui.ImageByteFormat.png);

    if (pngBytes != null) {
      final directory = await getApplicationDocumentsDirectory();
      final path = '${directory.path}/noise.png';
      final file = File(path);
      await file.writeAsBytes(pngBytes.buffer.asUint8List());
    }
  }

  Future<ui.Image> generateNoise() async {
    final ui.PictureRecorder recorder = ui.PictureRecorder();
    final Canvas canvas = Canvas(recorder);
    final paint = Paint();

    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        // double value = noise3_ImproveXY(seed, x * frequency, y * frequency, 0.0);
        double value = OpenSimplex2S()
            .noise3ImproveXY(seed, x * frequency, y * frequency, 0.0);
        int grayscale = ((value + 1) * 127.5).toInt();
        print(value);
        int rgb = 0xFF000000 | (grayscale << 16) | (grayscale << 8) | grayscale;
        paint.color = Color(rgb);
        canvas.drawRect(Rect.fromLTWH(x.toDouble(), y.toDouble(), 1, 1), paint);
      }
    }

    final ui.Picture picture = recorder.endRecording();
    final ui.Image image = await picture.toImage(width, height);
    return image;
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: Text('OpenSimplex Noise Demo'),
      ),
      body: Center(
        child: FutureBuilder<ui.Image>(
          future: generateNoise(),
          builder: (BuildContext context, AsyncSnapshot<ui.Image> snapshot) {
            if (snapshot.connectionState == ConnectionState.done &&
                snapshot.hasData) {
              return CustomPaint(
                painter: ImagePainter(snapshot.data!),
              );
            } else if (snapshot.hasError) {
              return Text('Error: ${snapshot.error}');
            } else {
              return CircularProgressIndicator();
            }
          },
        ),
      ),
      floatingActionButton: FloatingActionButton(
        onPressed: generateNoiseImage,
        child: Icon(Icons.save),
      ),
    );
  }
}

class ImagePainter extends CustomPainter {
  final ui.Image image;

  ImagePainter(this.image);

  @override
  void paint(Canvas canvas, Size size) {
    canvas.drawImage(image, Offset.zero, Paint());
  }

  @override
  bool shouldRepaint(CustomPainter oldDelegate) {
    return false;
  }
}

void main() {
  runApp(MaterialApp(
    home: OpenSimplexDemo(),
  ));
}

// Placeholder for the noise3_ImproveXY function
double noise3_ImproveXY(int seed, double x, double y, double z) {
  // Proper implementation of OpenSimplex2S noise algorithm needed here
  // For now, we use a dummy implementation
  return (Random(seed).nextDouble() * 2 - 1); // Dummy return value for now
}
