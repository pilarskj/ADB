package test.adbp;

import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

public class StreamTest {

    @Test
    public void testStreams() {

        int m = 1000;
        double t_or = 20;
        double dx = t_or / (m - 1);

        // loop
        /*
        double[] t_seq = new double[m];
        for (int i = 0; i < m; i++) {
            t_seq[i] = dx * i;
        } // time: 11 ms
        */

        // Stream
        // double[] t_seq = DoubleStream.iterate(0, n -> n + dx).limit(m).toArray();
        double[] t_seq = IntStream.range(0, m)
                .mapToDouble(i -> i * dx)
                .toArray(); // 13 ms

        // System.out.println(Arrays.toString(t_seq));
    }
}
