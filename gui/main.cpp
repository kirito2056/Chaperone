#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <QQmlContext>
#include <QUrl>

int main(int argc, char *argv[])
{
    QGuiApplication app(argc, argv);

    QQmlApplicationEngine engine;
    engine.rootContext()->setContextProperty(
        QStringLiteral("defaultPdbPath"), QStringLiteral(CHAPERONE_PDB_PATH));
    engine.rootContext()->setContextProperty(
        QStringLiteral("runsDir"), QStringLiteral(CHAPERONE_RUNS_DIR));

    QObject::connect(
        &engine, &QQmlApplicationEngine::objectCreationFailed, &app,
        []() { QCoreApplication::exit(1); }, Qt::QueuedConnection);

    engine.load(QUrl::fromLocalFile(QStringLiteral(CHAPERONE_QML_DIR "/Main.qml")));
    return app.exec();
}
