package main;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public class FASTQ {
    private static final int HEADER_LINE = 0;
    private static final int SEQUENCE_LINE = 1;
    private static final int PLUS_LINE = 2;
    private static final int QUALITY_LINE = 3;

    private final Path path;

    private long readingsCount = 0;
    private long totalReadingsLength = 0;
    private long gcCount = 0;
    private long minReadingLength = Long.MAX_VALUE;
    private long maxReadingLength = Long.MIN_VALUE;

    public FASTQ(String path) throws FASTQException {
        this.path = Path.of(path);

        computeReadingStatistic();
    }

    public long readingsCount() throws FASTQException {
        return readingsCount;
    }

    public long minReadingLength() throws FASTQException {
        return minReadingLength;
    }

    public long maxReadingLength() throws FASTQException {
        return maxReadingLength;
    }

    public double averageReadingLength() throws FASTQException {
        return (double) totalReadingsLength / readingsCount;
    }

    public double gcComposition() throws FASTQException {
        return (double) gcCount / totalReadingsLength * 100;
    }

    /**
     * Calculates the average quality value on the Phred scale at the specified position
     */
    public double phredAt(int position) throws FASTQException {
        int lineType = 0;
        int count = 0;
        int phredSum = 0;

        try (BufferedReader reader = Files.newBufferedReader(path)) {

            String line;
            while ((line = reader.readLine()) != null) {
                if (lineType == QUALITY_LINE && line.length() >= position) {
                    count++;
                    char ch = line.charAt(position - 1);
                    phredSum += (int) ch - 33;
                }

                lineType = ((lineType + 1) % 4);
            }
        } catch (IOException e) {
            throw new FASTQException(e);
        }

        return count > 0 ? (double) phredSum / count : 0;
    }

    public long trim(int windowLength, int targetQuality) throws FASTQException {
        Path outputPath = generateOutputPath(windowLength, targetQuality);
        long trimmedCount = 0; // count of deleted sequences after trimming

        try (BufferedReader reader = Files.newBufferedReader(path);
             BufferedWriter writer = Files.newBufferedWriter(outputPath)) {

            String[] read;
            while ((read = readNext(reader)) != null) {
                String sequenceLine = read[SEQUENCE_LINE];
                String qualityLine = read[QUALITY_LINE];

                String[] trimmed = trim(sequenceLine, qualityLine, windowLength, targetQuality); // [seq, qual]

                if (trimmed == null || trimmed[0].isEmpty()) {
                    trimmedCount++;
                    continue;
                }

                if (trimmed[0].length() < sequenceLine.length()) {
                    String updatedPlusLine = updateLengthInPlusLine(read[PLUS_LINE], trimmed[0].length());
                    String[] toWrite = new String[]{read[HEADER_LINE], trimmed[0], updatedPlusLine, trimmed[1]};

                    write(writer, toWrite);
                } else {
                    write(writer, read);
                }
            }
        } catch (IOException e) {
            throw new FASTQException(e);
        }

        return trimmedCount;
    }

    /**
     * Updates value of length in PLUS_LINE
     */
    private String updateLengthInPlusLine(String plusLine, int newLength) {
        if (plusLine.contains("length=")) {
            return plusLine.replaceAll("length=\\d+", "length=" + newLength);
        }

        return plusLine;
    }


    /**
     * Write 4 lines describing one reading
     */
    private void write(BufferedWriter writer, String[] reading) throws IOException {
        if (reading.length != 4) return;

        for (String s : reading) {
            writer.write(s);
            writer.newLine();
        }
    }

    private Path generateOutputPath(int windowLength, int quality) {
        String fileName = path.getFileName().toString();
        String newFileName = String.format("%s.tr_w%d_q%d", fileName, windowLength, quality);

        return path.resolveSibling(newFileName);
    }

    private String[] trim(String sequenceLine, String qualityLine, int windowLength, int requiredQuality) {
        int totalRequiredQuality = requiredQuality * windowLength;

        int[] quals = getQualityAsInteger(qualityLine);

        if (quals.length < windowLength) {
            return null;
        }

        int total = 0;
        for (int i = 0; i < windowLength; i++)
            total += quals[i];

        if (total < totalRequiredQuality)
            return null;

        int lengthToKeep = quals.length;

        for (int i = 0; i < quals.length - windowLength; i++) {
            total = total - quals[i] + quals[i + windowLength];
            if (total < totalRequiredQuality) {
                lengthToKeep = i + windowLength;
                break;
            }
        }

        int i = lengthToKeep;

        if (i < 1)
            return null;

        if (i < quals.length) {
            return new String[]{
                    sequenceLine.substring(0, i),
                    qualityLine.substring(0, i)
            };
        }

        return new String[]{sequenceLine, qualityLine};
    }

    private int[] getQualityAsInteger(String qualityLine) {
        int[] qual = new int[qualityLine.length()];

        for (int i = 0; i < qualityLine.length(); i++) {
            char ch = qualityLine.charAt(i);
            qual[i] = (int) ch - 33;
        }

        return qual;
    }

    /**
     * Returns 4 lines describing next sequence
     */
    private String[] readNext(BufferedReader reader) throws IOException {
        String[] read = new String[4];

        if ((read[0] = reader.readLine()) != null) {
            for (int i = 1; i < read.length; i++) {
                read[i] = reader.readLine();
            }
            return read;
        }

        return null;
    }

    private void computeReadingStatistic() throws FASTQException {
        int lineType = 0;

        try (BufferedReader reader = Files.newBufferedReader(path)) {
            String line;

            while ((line = reader.readLine()) != null) {
                if (lineType == HEADER_LINE) {
                    this.readingsCount++;
                }

                if (lineType == SEQUENCE_LINE) {
                    for (int i = 0; i < line.length(); i++) {
                        char ch = line.charAt(i);
                        if (ch == 'G' || ch == 'C') this.gcCount++;
                    }

                    int length = line.length();
                    this.totalReadingsLength += length;
                    this.minReadingLength = Math.min(this.minReadingLength, length);
                    this.maxReadingLength = Math.max(this.maxReadingLength, length);
                }

                lineType = ((lineType + 1) % 4);
            }
        } catch (IOException e) {
            throw new FASTQException(e);
        }
    }

    public static class FASTQException extends Exception {
        public FASTQException(Exception cause) {
            super(cause);
        }
    }
}
