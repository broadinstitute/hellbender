package org.broadinstitute.hellbender.utils.io;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;

/**
 * Class to hold a set of {@link Path} to be delete on the JVM exit through a shutdown hook.
 *
 * <p>This class is a modification of {@link htsjdk.samtools.util.nio}
 */
class DeleteRecursivelyOnExitPathHook {
    private static LinkedHashSet<Path> paths = new LinkedHashSet<>();
    static {
        Runtime.getRuntime().addShutdownHook(new Thread(DeleteRecursivelyOnExitPathHook::runHooks));
    }

    private DeleteRecursivelyOnExitPathHook() {}

    /**
     * Adds a {@link Path} for deletion on JVM exit.
     *
     * @param path path to be deleted.
     *
     * @throws IllegalStateException if the shutdown hook is in progress.
     */
    public static synchronized void add(Path path) {
        if(paths == null) {
            // DeleteOnExitHook is running. Too late to add a file
            throw new IllegalStateException("Shutdown in progress");
        }

        paths.add(path);
    }

    static void runHooks() {
        LinkedHashSet<Path> thePaths;

        synchronized (DeleteRecursivelyOnExitPathHook.class) {
            thePaths = paths;
            paths = null;
        }

        ArrayList<Path> toBeDeleted = new ArrayList<>(thePaths);

        // reverse the list to maintain previous jdk deletion order.
        // Last in first deleted.
        Collections.reverse(toBeDeleted);
        for (Path path : toBeDeleted) {
            try {
                IOUtils.deleteRecursively(path);
            } catch (IOException | SecurityException e) {
                // do nothing if cannot be deleted, because it is a shutdown hook
            }
        }
    }
}