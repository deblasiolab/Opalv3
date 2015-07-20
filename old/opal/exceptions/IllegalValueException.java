package opal.exceptions;

/**
 * Exception class for access in empty containers
 * such as stacks, queues, and priority queues.
 * @author Mark Allen Weiss
 */
public class IllegalValueException extends RuntimeException {
    /**
     * Construct this exception object.
     * @param message the error message.
     */
    public IllegalValueException( String message ) {
        super( message );
    }
}