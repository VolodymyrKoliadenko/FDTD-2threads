
import java.awt.*;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import javax.swing.JFrame;

public class FirstWindow extends JFrame {

    Label ny, nx, widthL;
    public Label stroka;
    Button start;

    public FirstWindow() {

        super("First FDTD window");
        setSize(500, 650);
        setResizable(false);

        TextField textNy = new TextField("500", 1);
        textNy.setBounds(0, 0, 40, 21);
        textNy.setBackground(Color.YELLOW);
        textNy.addTextListener(
                (te -> {
                    String s = textNy.getText();
                    try {
                        double t = Double.parseDouble(s);
                        int t2 = (int) Math.round(t);;
                        BasicEx.ny = t2;
                        textNy.setBackground(new Color(200, 255, 200));
                        ny.setText("Ny = " + t2 + " pixels");
                    } catch (Exception e) {
                        textNy.setBackground(new Color(255, 100, 100));
                    }
                })
        );
        add(textNy);
        ny = new Label("Ny");
        ny.setBounds(41, 0, 190, 21);
        add(ny);

        TextField textNx = new TextField("1024", 1);
        textNx.setBounds(0, 30, 40, 21);
        textNx.setBackground(Color.YELLOW);
        textNx.addTextListener(
                (te -> {
                    String s = textNx.getText();
                    try {
                        double t = Double.parseDouble(s);
                        int t2 = (int) round(log(t) / log(2));
                        int t3 = (int) pow(2, t2);
                        BasicEx.nx = t3;
                        textNx.setBackground(new Color(200, 255, 200));
                        nx.setText("Nx = " + t3 + " pixels");
                    } catch (Exception e) {
                        textNx.setBackground(new Color(255, 100, 100));
                    }
                })
        );
        add(textNx);
        nx = new Label("Nx");
        nx.setBounds(41, 30, 190, 21);
        add(nx);

        TextField textQ = new TextField("1.08", 1);
        textQ.setBounds(0, 60, 40, 21);
        textQ.setBackground(Color.YELLOW);
        textQ.addTextListener(
                (te -> {
                    String s = textQ.getText();
                    try {
                        double t = Double.parseDouble(s);
                        FDTD.Q = (float) t;
                        textQ.setBackground(new Color(200, 255, 200));
                    } catch (Exception e) {
                        textQ.setBackground(new Color(255, 100, 100));
                    }
                })
        );
        add(textQ);
        Label q = new Label("Q");
        q.setBounds(42, 60, 190, 21);
        add(q);

        TextField textM = new TextField("0.016", 1);
        textM.setBounds(0, 90, 40, 21);
        textM.setBackground(Color.YELLOW);
        textM.addTextListener(
                (te -> {
                    String s = textM.getText();
                    try {
                        double t = Double.parseDouble(s);
                        FDTD.mod = (float) t;
                        textM.setBackground(new Color(200, 255, 200));
                    } catch (Exception e) {
                        textM.setBackground(new Color(255, 100, 100));
                    }
                })
        );
        add(textM);
        Label mod = new Label("Δε");
        mod.setBounds(42, 90, 190, 21);
        add(mod);

        TextField textW = new TextField("1.9", 1);
        textW.setBounds(0, 120, 40, 21);
        textW.setBackground(Color.YELLOW);
        textW.addTextListener(
                (te -> {
                    String s = textW.getText();
                    try {
                        double t = Double.parseDouble(s);
                        FDTD.w = (float) (t / 1e6);
                        textW.setBackground(new Color(200, 255, 200));
                    } catch (Exception e) {
                        textW.setBackground(new Color(255, 100, 100));
                    }
                })
        );
        add(textW);
        Label w = new Label("Gauss half-width, µm");
        w.setBounds(42, 120, 190, 21);
        add(w);

        TextField textLambda = new TextField("1064", 1);
        textLambda.setBounds(0, 150, 40, 21);
        textLambda.setBackground(Color.YELLOW);
        textLambda.addTextListener(
                (te -> {
                    String s = textLambda.getText();
                    try {
                        double t = Double.parseDouble(s);
                        FDTD.lambd = (float) (t / 1e9);
                        textLambda.setBackground(new Color(200, 255, 200));
                    } catch (Exception e) {
                        textLambda.setBackground(new Color(255, 100, 100));
                    }
                })
        );
        add(textLambda);
        Label lambda = new Label("Wawelength, nm");
        lambda.setBounds(42, 150, 190, 21);
        add(lambda);

        TextField textWidth = new TextField("73", 1);
        textWidth.setBounds(0, 180, 40, 21);
        textWidth.setBackground(Color.YELLOW);
        textWidth.addTextListener(
                (te -> {
                    String s = textWidth.getText();
                    try {
                        double t0 = Double.parseDouble(s);
                        double t = t0 / 1e6 / BasicEx.nx;
                        FDTD.dx = (float) t;
                        double deLambd = round(FDTD.lambd / t * 100) / 100.0;
                        widthL.setText("Width " + (int) t0 + " µm. "
                                + deLambd + " px/λ");
                        if (deLambd < 10) {
                            widthL.setBackground(new Color(255, 100, 100));
                        } else if (deLambd < 15) {
                            widthL.setBackground(new Color(255, 255, 100));
                        } else {
                            widthL.setBackground(new Color(100, 255, 100));
                        }
                        textWidth.setBackground(new Color(200, 255, 200));
                    } catch (Exception e) {
                        textWidth.setBackground(new Color(255, 100, 100));
                    }
                })
        );
        add(textWidth);
        widthL = new Label("Width, µm");
        widthL.setBounds(42, 180, 150, 21);
        add(widthL);

        start = new Button("Start calculating");  // ne pashut??
        start.addActionListener(
                (ae -> {
                    //this.remove(start);
                    start.setLabel("Next calculation");
                    //remove(start);
                    BasicEx ex = new BasicEx();
                    ex.setVisible(true);

                })
        );
        start.setBounds(0, 400, 100, 40); //BUTTON
        add(start);

        stroka = new Label("", Label.CENTER);
        stroka.setFont(new Font("Arial", Font.BOLD, 20));
        stroka.setBounds(109, 592, 299, 23);
        add(stroka);

        //button exit !
        setLayout(null);
        setVisible(true);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

//        FileDialog fDialog = new FileDialog(frame, "Choose image path", FileDialog.SAVE);
//        fDialog.setVisible(true);
//        while (fDialog.getFile() == null) {
//        };
//        File imageFile = new File(fDialog.getDirectory(), fDialog.getFile());

    }

//    private JFrame thiss() {
//        return this;
//    }
}
