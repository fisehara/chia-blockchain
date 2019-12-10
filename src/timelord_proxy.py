import asyncio
import logging

log = logging.getLogger(__name__)


async def redirect_stream(reader, writer):
    while True:
        try:
            chunk = await reader.read(8)
            if (chunk == b''):
                print("Empty buffer detected... stopping")
                break
            writer.write(chunk)
            await writer.drain()
        except Exception as e:
            print(f"Caught exception: {e}... stopping")
            break


async def handle_echo(reader, writer):
    ip = '127.0.0.1'
    port = await reader.readexactly(4)
    str_port = port.decode()

    proc = await asyncio.create_subprocess_shell(
        f"./lib/chiavdf/fast_vdf/vdf_server {str_port}"
    )

    # Communication between the proxy and the spawned process.
    proxy_reader, proxy_writer = None, None

    for _ in range(10):
        try:
            proxy_reader, proxy_writer = await asyncio.open_connection(ip, str_port)
            socket = proxy_writer.get_extra_info("socket")
            socket.settimeout(None)
            break
        except Exception as e:
            log.info(f"Exception while connecting to the process: {str(e)}")
            await asyncio.sleep(1)

    if not proxy_reader or not proxy_writer:
        raise Exception("Unable to connect to VDF server")

    event_loop = asyncio.get_event_loop()
    t1 = event_loop.create_task(redirect_stream(reader, proxy_writer))
    t2 = event_loop.create_task(redirect_stream(proxy_reader, writer))

    await proc.wait()
    t1.cancel()
    t2.cancel()

loop = asyncio.get_event_loop()
coro = asyncio.start_server(handle_echo, '127.0.0.1', 9999, loop=loop)
server = loop.run_until_complete(coro)

# Serve requests until Ctrl+C is pressed
print('Serving on {}'.format(server.sockets[0].getsockname()))
try:
    loop.run_forever()
except KeyboardInterrupt:
    pass

# Close the server
server.close()
loop.run_until_complete(server.wait_closed())
loop.close()
