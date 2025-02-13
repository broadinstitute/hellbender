package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.concurrent.LinkedBlockingQueue;

public abstract class ThreadedReadWalker extends ReadWalker implements Runnable {

    public static final int CAPACITY = 5000;

    // apply() message and queue - carrying read/referenceContext tuples
    static class ApplyMessage {
        GATKRead read;
        ReferenceContext referenceContext;
        FeatureContext featureContext;
    }

    private LinkedBlockingQueue<ApplyMessage> applyQueue;
    private Thread  worker;

    /**
     * turn threaded walker on
     **/
    @Argument(fullName = "threaded-walker", doc = "turn threaded walker on?", optional = true)
    public boolean threadedWalker = false;


    @Override
    public void onTraversalStart() {

        // if running in threaded mode, start worker thread
        if ( threadedWalker ) {
            applyQueue = new LinkedBlockingQueue<>(CAPACITY);
            worker = new Thread(this);
            worker.setUncaughtExceptionHandler(getUncaughtExceptionHandler());
            worker.start();
        }
    }

    @Override
    public void closeTool() {

        // if running in threaded mode, signal end of input and wait for thread to complete
        if ( threadedWalker ) {
            try {
                applyQueue.put(new ApplyMessage());
                worker.join();
            } catch (InterruptedException e) {
                throw new GATKException("", e );
            }
        }
    }

    @Override
    final public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        if ( !acceptRead(read, referenceContext, featureContext) ) {
            return;
        }

        if ( threadedWalker ) {
            // redirect message into queue
            ApplyMessage message = new ApplyMessage();
            message.read = read;
            message.referenceContext = referenceContext;
            message.featureContext = featureContext;
            try {
                applyQueue.put(message);
            } catch (InterruptedException e) {
                throw new GATKException("", e);
            }
        } else {
            applyPossiblyThreaded(read, referenceContext, featureContext);
        }
    }

    @Override
    public void run() {

        // loop on messages
        ApplyMessage message;
        while ( true ) {
            // take next message off the queue
            try {
                message = applyQueue.take();
            } catch (InterruptedException e) {
                throw new GATKException("", e);
            }

            // end of stream?
            if ( message.read == null ) {
                break;
            }

            // process message
            applyPossiblyThreaded(message.read, message.referenceContext, message.featureContext);
        }
    }

    public abstract boolean acceptRead(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext);
    public abstract void applyPossiblyThreaded(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext);

    public Thread.UncaughtExceptionHandler getUncaughtExceptionHandler() {
        Thread.UncaughtExceptionHandler handler = (th, ex) -> {
            System.out.println("Uncaught exception: " + ex);
            ex.printStackTrace(System.out);
            System.exit(-1);
        };
        return handler;
    }
}
