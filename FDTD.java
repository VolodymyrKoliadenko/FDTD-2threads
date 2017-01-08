
import static java.lang.Math.*;
//import java.util.Calendar;

public class FDTD {

    public static float lambd = 1064e-9f;
    public static float dx = lambd / 15;         // max dx=lam/10 !!
    public static float Q = 1.08f;//1.0 - center
    public static float mod = 0.008f * 2;//max modulation of EPSILON =2*deltaN;
    public static float w = 19e-7f;

    public static float[][] callFDTD(int nx, int ny, String method) {
        int i, j;//for loops
        float x; //for coordinate
        final float dy = dx;
        final float period = 2e-6f;
        final float n = 1;//refractive index

        final float prodol = 2 * n * period * period / lambd / Q;

        float[][] Ez = new float[nx + 1][ny + 1];
        float[][] Hx = new float[nx - 1 + 1][ny - 1 + 1];
        float[][] Hy = new float[nx][ny];

        final float omega = (float) (2 * PI / lambd);
        final float dt = dx / 2;
        final float tau = 2e-5f * 999;// light decay
        //final float Z = 376.7303f;
        final float s = dt / dx;
        final float k3 = (1 - s) / (1 + s);//for MUR
        final float alpha = (float) (sin(0.0 / 180 * PI));//radians

        final int begin = 10;// must bee small
        System.out.println("Imeem__Sloy= " + (ny - begin) * dy / prodol * 2);

        final float ds = dt * dt / dx / dx;// renormalization  dH=dt/dx/Z;
        float[][] e = new float[nx + 1][ny + 1];
        for (i = 1; i < nx + 1; i++) {
            for (j = 1; j < ny + 1; j++) {
                e[i][j] = (float) (ds / (n + ((j < begin) ? 0 : (mod / 2) * (1 + signum(-0.1 
                        + cos(2 * PI * (i - nx / 2.0 + 0.5) * dx / period) * sin(2 * PI * (j - begin) * dy / prodol))))));
            }
        }

        float[][] end = new float[2][nx + 1]; // boundary conditions
        float[][] top = new float[2][ny + 1];
        float[][] bottom = new float[2][ny + 1];
        
        final int tMax = (int) (ny * 2.2);
        System.out.println("START CICLE");
        for (int t = 1; t <= tMax; t++) {                    // begin main loop
            if ((t % 10 == 0) && (method.equals("cos"))) {
                BasicEx.fW.stroka.setText("Calculated " + t + " of " + tMax + " steps");
            }
            float tt = Math.min(t * s + 10, ny - 1);
            //gauss
            switch (method) {
                case "cos":        
                    for (i = 1; i <= nx - 1; i++) {
                        x = (float) (dx * (i - (float) nx / 2 + 0.5));
                        Ez[i][1] = (float) (exp(-pow(x, 2) / w / w - (t - 1) * dt / tau)
                                * cos((x * alpha + (t - 1) * dt) * omega));
                    }
                    break;
                case "sin":
                    for (i = 1; i <= nx - 1; i++) {
                        x = (float) (dx * (i - (float) nx / 2 + 0.5));
                        Ez[i][1] = (float) (exp(-pow(x, 2) / w / w - (t - 1) * dt / tau) * sin((x * alpha + (t - 1) * dt) * omega));
                    }
                    break;
            }

            for (i = 1; i <= nx; i++) {  // boundary conditions
                end[0][i] = Ez[i][ny - 1];
                end[1][i] = Ez[i][ny];
            }
            System.arraycopy(Ez[1], 0, top[0], 0, ny + 1);
            System.arraycopy(Ez[2], 0, top[1], 0, ny + 1);
            System.arraycopy(Ez[nx - 1], 0, bottom[0], 0, ny + 1);
            System.arraycopy(Ez[nx], 0, bottom[1], 0, ny + 1);

            for (i = 2; i <= nx - 1; i++) {        // main Ez  
                for (j = 2; j <= tt; j++) {   
                    Ez[i][j] += e[i][j] * ((Hx[i][j - 1] - Hx[i][j] + Hy[i][j] - Hy[i - 1][j]));
                }
            }

            for (i = 1; i <= nx; i++) {    // boundary conditions
                Ez[i][ny] = end[0][i] + k3 * (end[1][i] - Ez[i][ny - 1]);//end
            }
            for (i = 1; i <= ny; i++) {
                Ez[1][i] = top[1][i] + k3 * (top[0][i] - Ez[2][i]);
                Ez[nx][i] = bottom[0][i] + k3 * (bottom[1][i] - Ez[nx - 1][i]);
            }
           switch (method) {
                case "cos":
                    for (i = 1; i <= nx - 1; i++) {
                        x = (float) (dx * (i - (float) nx / 2 + 0.5));
                        Ez[i][1] = (float) (exp(-pow(x, 2) / w / w - (t - 1) * dt / tau) * cos((x * alpha + t * dt) * omega));
                    }
                    break;
                case "sin":
                    for (i = 1; i <= nx - 1; i++) {
                        x = (float) (dx * (i - (float) nx / 2 + 0.5));
                        Ez[i][1] = (float) (exp(-pow(x, 2) / w / w - (t - 1) * dt / tau) * sin((x * alpha + t * dt) * omega));
                    }
                    break;
            }
           
            for (i = 1; i <= nx - 1; i++) {        // main Hx Hy
                for (j = 1; j <= tt; j++) {
                    Hx[i][j] += Ez[i][j] - Ez[i][j + 1];
                    Hy[i][j] += Ez[i + 1][j] - Ez[i][j];
                }
            }
        }

        int pos = method.equals("cos") ? 0 : 1;
        BasicEx.forFurier[pos] = new double[nx];
        int endF = (int) (ny * 0.95);//0.99
        for (i = 1; i <= nx; i++) {
            BasicEx.forFurier[pos][i - 1] = Ez[i][endF];
            for (j = 1; j <= ny; j++) {
                Ez[i][j] = abs(Ez[i][j]);// ABS
            }
        }
        Hx = null;
        Hy = null;
        e = null;
        return Ez;

    }

}
