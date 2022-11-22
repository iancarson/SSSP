/*
    Coordinator.java

    Serves to slow down execution of a collection of threads, so that
    behavior is visible on the screen, and to notify all running threads
    when the user wants them to die.

    (c) Michael L. Scott, 2022; based on code originally written in 2002.
    For use by students in CSC 2/454 at the University of Rochester,
    during the Fall 2022 term.  All other use requires written permission
    of the author.
 */

// import java.util.*;
// import java.lang.*;

class Coordinator {
    private boolean open = true;
    // set to false temporarily when threads are supposed to die.
    private boolean running = true;
    // set to false when threads are supposed to pause.
    private int numThreads = 0;
    // number of active worker threads.  Maintained by register and
    // unregister methods.

    // A thread terminates early by throwing itself a KilledException.
    //
    public class KilledException extends Throwable {}

    public synchronized void register() {
        numThreads++;
    }

    public synchronized void unregister() {
        numThreads--;
        notifyAll();
        // so event thread knows to inspect numThreads again
    }

    // Pause or die if so instructed.
    //
    private synchronized void gate()
            throws KilledException {
        if (!open) {
            numThreads--;
            notify();
            throw new KilledException();
        }
        while (!running) {
            try {
                wait();
            } catch(InterruptedException e) {};
            if (!open) {
                numThreads--;
                notify();
                throw new KilledException();
            }
        }
    }

    // Wait a bit before proceeding through gate.
    //
    public void hesitate() throws KilledException {
        try {
            Thread.sleep(50);   // milliseconds
        } catch(InterruptedException e) {};
        gate();
    }

    // Toggle running.  Resume paused threads if appropriate.
    //
    public synchronized void toggle() {
        running = !running;
        if (running) {
            notifyAll();
        }
    }

    // Kill all threads using the coordinator.
    //
    public synchronized void stop() {
        open = false;
        notifyAll();
        while (numThreads > 0) {
            try {
                wait();
            } catch(InterruptedException e) {};
        }
        open = true;
        running = true;
    }

    // (Default constructor only)
}
