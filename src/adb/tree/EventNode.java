package adb.tree;

import beast.base.core.Description;

@Description("A node marking a branching event along the lineage.")
public class EventNode {

    protected int type;
    protected double height;

    public EventNode(int type, double height) {
        this.type = type;
        this.height = height;
    }

    public int getType() {
        return type;
    }

    public void setType(int type) {
        this.type = type;
    }

    public double getHeight() {
        return height;
    }

    public void setHeight(double height) {
        this.height = height;
    }

}
